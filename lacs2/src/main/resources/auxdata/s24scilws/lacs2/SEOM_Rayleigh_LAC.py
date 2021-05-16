# coding: utf-8
import os, sys
import gdal
import numpy as np
import scipy.ndimage
from PIL import Image
import datetime
import numpy.ma as ma

import subprocess
import os.path
import re
import time
from pysolar.solar import *
import copy
#from scipy.interpolate import Rbf
import pygrib

import glob
import scipy.interpolate as interpolate
import shutil
import pytz
from math import *
from math import fabs

import logging

logging.basicConfig()
logger = logging.getLogger('logger')


#Module pour interpoler un tableau en fonction de la distance aux points
def InterpolArray(array):
    arrayout=copy.deepcopy(array)
    for i in range(0,len(array)):
        for j in range(0,len(array[0])):
            if np.isnan(array[i,j])==True:
                NbPoint=0
                SommePoids=0
                SommeVal=0
                Delta=1
                while NbPoint<4 and Delta<max(len(array),len(array[0])):
                    for ki in range(max(0,i-Delta),min(len(array),i+Delta+1)):
                        for kj in range(max(0,j-Delta),min(len(array[0]),j+Delta+1)):
                            if (abs(ki-i)==Delta or abs(kj-j)==Delta) and np.isnan(array[ki,kj])==False:
                                Distance=1.0*((ki-i)**2+(kj-j)**2)
                                SommePoids=SommePoids+1/Distance
                                SommeVal=SommeVal+array[ki,kj]/Distance
                                NbPoint=NbPoint+1
                    Delta=Delta+1
                arrayout[i,j]=SommeVal/SommePoids
    return arrayout

#Module pour les concentrations en O3
def getO3(Dir):
    ListfileECMWFT=glob.glob(os.path.join(Dir, 'AUX_DATA', 'AUX_ECMWFT'))


    if len(ListfileECMWFT)>0:
        fileECMWFT=ListfileECMWFT[0]
        datum = Dir.split('_')[-1].split('T')[0]

        grbs =pygrib.open("%s"%fileECMWFT)



        tvar = grbs.select(name="Total column ozone")[0]
        data, lats, lons = tvar.data(lat1=-100, lat2=100, lon1=-400, lon2=400)
        mean=data.reshape(1,len(data)*len(data[0]),1).mean()
        logging.debug ("{} {}".format(datum,mean))


        if os.path.isfile(os.path.join(newfolder, "MeanO3.txt")):
            os.remove(os.path.join(newfolder, "MeanO3.txt"))

        f = open(os.path.join(newfolder, 'MeanO3.txt'),'a+')
        f.write("%s;%s;" %(datum,mean))
        f.close


        fileB1=glob.glob(os.path.join(Dir, 'IMG_DATA', '*B01.jp2'))

        inDs = gdal.Open("%s" %fileB1[0])
        rows = inDs.RasterYSize
        cols = inDs.RasterXSize
        myband = inDs.GetRasterBand(1)
        driver = gdal.GetDriverByName('GTiff')
        inGs=inDs.GetGeoTransform()
        GeoT=[inGs[0],rows*inGs[1]/len(data),inGs[2],inGs[3],inGs[4],cols*inGs[5]/len(data[0])]


        outDs = driver.Create(os.path.join(newfolder, f"O3-{datum}.TIF"), len(data[1]), len(data), 1, gdal.GDT_Float32, options=['TILED=YES','BIGTIFF=IF_SAFER','BLOCKXSIZE=512','BLOCKYSIZE=512','COMPRESS=DEFLATE'])
        outBand = outDs.GetRasterBand(1)
        outBand.WriteArray(data)
        outBand.FlushCache()
        outBand.SetNoDataValue(-99)
        outDs.SetGeoTransform(GeoT)
        outDs.SetProjection(inDs.GetProjection())

        mean=mean*46698 #Conversion en Dobson
        grbs.close()
    else:
        mean=332.8

    return mean

#Module pour les pressions
def getPressure(Dir):
    ListfileECMWFT=glob.glob(os.path.join(Dir, 'AUX_DATA', 'ECMWFT'))
    ArrayPres=[]

    if len(ListfileECMWFT)>0:
        fileECMWFT=ListfileECMWFT[0]
        tmp=fileECMWFT.split("/")
        datum=tmp[len(tmp)-3]
        grbs =pygrib.open("%s" %fileECMWFT)

        tvar = grbs.select(name="Mean sea level pressure")[0]
        ArrayPres, lats, lons = tvar.data(lat1=-100, lat2=100, lon1=-400, lon2=400)
        mean=ArrayPres.reshape(1,len(ArrayPres)*len(ArrayPres[0]),1).mean()
        logging.debug ("{} {}".format(datum,mean))


        fileB1=glob.glob(os.path.join(Dir, 'IMG_DATA', '*B01.jp2'))

        inDs = gdal.Open("%s" %fileB1[0])
        rows = inDs.RasterYSize
        cols = inDs.RasterXSize
        myband = inDs.GetRasterBand(1)
        driver = gdal.GetDriverByName('GTiff')
        inGs=inDs.GetGeoTransform()
        GeoT=[inGs[0],rows*inGs[1]/len(ArrayPres),inGs[2],inGs[3],inGs[4],cols*inGs[5]/len(ArrayPres[0])]


        outDs = driver.Create(os.path.join(Dir, f"Pressure-{datum}.TIF"), len(ArrayPres[1]), len(ArrayPres), 1, gdal.GDT_Float32, options=['TILED=YES','BIGTIFF=IF_SAFER','BLOCKXSIZE=512','BLOCKYSIZE=512','COMPRESS=DEFLATE'])
        outBand = outDs.GetRasterBand(1)
        outBand.WriteArray(ArrayPres)
        outBand.FlushCache()
        outBand.SetNoDataValue(-99)
        outDs.SetGeoTransform(GeoT)
        outDs.SetProjection(inDs.GetProjection())

        grbs.close()

    return ArrayPres

#Module pour recuperer absorbance O3
def getTau(Ly):
    Lambda=float(Ly)*1000

    ArrayTau=[[0,0],[412.5,0.0002179],[442.5,0.002814],[490,0.02006],[510,0.04081],[560,0.104],[620,0.109],[665,0.0505],[681.25,0.03526],[708.75,0.01881],[753.75,0.008897],[761.875,0.006634],[778.75,0.007693],[865,0.002192],[885,0.001211],[900,0.001517],[1200,0],[20000,0]]

    i=0
    while i < len(ArrayTau):

        if float(ArrayTau[i][0])<Lambda:
            LambdaMoins=float(ArrayTau[i][0])
            TauMoins=float(ArrayTau[i][1])
            LambdaPlus=float(ArrayTau[i+1][0])
            TauPlus=float(ArrayTau[i+1][1])
            Tau=(Lambda-LambdaMoins)*(TauPlus-TauMoins)/(LambdaPlus-LambdaMoins)+TauMoins

        i=i+1

    return Tau

#Recuperation des angles de visees du satellite
def getAngles2(filename):
    fic=open("%s" %filename,"r")
    lineF = fic.readline().strip()
    cpt=0
    bool=0
    ArrayViewingAngleZen=[]
    ArrayViewingAngleAzi=[]


    for i in range(0,20):
        tmp=[]
        tmp2=[]
        for j in range(0,20):
            tmp.append([])
            tmp2.append([])
        ArrayViewingAngleZen.append(tmp)
        ArrayViewingAngleAzi.append(tmp2)



    while "<n1:Quality_Indicators_Info" not in  lineF :
        lineF = fic.readline().strip()

        if "SENSING_TIME" in lineF:
            bool=bool+1
            TimeStamp=(((lineF.split(">"))[1]).split("<"))[0]
        if "<Viewing_Incidence_Angles_Grids bandId=\"2\"" in lineF:
            bool=bool+1
            lineF = fic.readline().strip()
            lineF = fic.readline().strip()
            lineF = fic.readline().strip()
            lineF = fic.readline().strip()
            lineF = fic.readline().strip()
            cpti=0
            for i in range(0,23):
                if i==5 or i ==11 or i==17:
                    lineF = fic.readline().strip()
                    continue

                tmp=lineF.split(">")
                tmp2=tmp[1].split("<")
                tmp3=tmp2[0].split(" ")
                cptj=0
                for j in range(0,23):
                    if j==5 or j==11 or j==17:
                        continue

                    if tmp3[j] != "NaN":
                        ArrayViewingAngleZen[cpti][cptj].append(float(tmp3[j]))
                    cptj=cptj+1
                lineF = fic.readline().strip()
                cpti=cpti+1

            lineF = fic.readline().strip()
            lineF = fic.readline().strip()
            lineF = fic.readline().strip()
            lineF = fic.readline().strip()
            lineF = fic.readline().strip()
            lineF = fic.readline().strip()

            cpti=0
            for i in range(0,23):
                if i==5 or i ==11 or i==17:
                    lineF = fic.readline().strip()
                    continue

                tmp=lineF.split(">")
                tmp2=tmp[1].split("<")
                tmp3=tmp2[0].split(" ")
                cptj=0
                for j in range(0,23):
                    if j==5 or j==11 or j==17:
                        continue
                    if tmp3[j] != "NaN":
                        ArrayViewingAngleAzi[cpti][cptj].append(float(tmp3[j]))
                    cptj=cptj+1
                lineF = fic.readline().strip()
                cpti=cpti+1

        if "Mean_Viewing_Incidence_Angle bandId=" in lineF:
            bool=bool+1
            lineF = fic.readline().strip()
            Zenith=(((lineF.split(">"))[1]).split("<"))[0]
            lineF = fic.readline().strip()
            Azimuth=(((lineF.split(">"))[1]).split("<"))[0]

        cpt=cpt+1
    fic.close()

    ViewingAngleZen=np.zeros((20,20))
    ViewingAngleAzi=np.zeros((20,20))

    for i in range(0,20):
        for j in range(0,20):
            ViewingAngleZen[i][j]=np.mean(ArrayViewingAngleZen[i][j])
            ViewingAngleAzi[i][j]=np.mean(ArrayViewingAngleAzi[i][j])
    '''
    x = np.random.rand(100)*4.0-2.0

    print ViewingAngleZen
    print "BBB"

    # 2-d tests - setup scattered data
    x = np.random.rand(100)*4.0-2.0
    y = np.random.rand(100)*4.0-2.0
    z = x*np.exp(-x**2-y**2)
    ti = np.linspace(-2.0, 2.0, 100)
    XI, YI = np.meshgrid(ti, ti)

    # use RBF
    rbf = Rbf(x, y, z, epsilon=2)
    ZI = rbf(XI, YI)
    print ZI


    ti = np.linspace(0, 19, 20)
    print ti
    XI, YI = np.meshgrid(ti, ti)

    print listX
    print listY
    print listZ
    rbf = Rbf(np.array(listX), np.array(listY), np.array(listZ),function='linear',epsilon=2)
    ZI = rbf(XI, YI)

    '''

    for i in range(0,20):
        if np.isnan(ViewingAngleZen[i][0])==False and np.isnan(ViewingAngleZen[i][19])==False:
            continue
        if np.isnan(ViewingAngleZen[i][0])==True:
            for j in range(0,20):
                if np.isnan(ViewingAngleZen[i][j])==False:
                    Slope=ViewingAngleZen[i][j]-ViewingAngleZen[i][j+1]
                    for k in range(0,j):
                        ViewingAngleZen[i][k]=ViewingAngleZen[i][j]+Slope*(j-k)
                    break
        if np.isnan(ViewingAngleZen[i][19])==True:
            for j in range(19,-1,-1):
                if np.isnan(ViewingAngleZen[i][j])==False:
                    Slope=ViewingAngleZen[i][j]-ViewingAngleZen[i][j-1]
                    for k in range(19,j,-1):
                        ViewingAngleZen[i][k]=ViewingAngleZen[i][j]-Slope*(j-k)
                    break

    ViewingAngleZen2=copy.deepcopy(ViewingAngleZen)
    for i in range(0,20):
            for j in range(0,20):
                if np.isnan(ViewingAngleZen[i][j])==True:
                    cpt=0
                    Somme=0
                    for ki in range(max(0,i-1),min(20,i+1)):
                        for kj in range(max(0,j-1),min(20,j+1)):
                            if np.isnan(ViewingAngleZen[ki][kj])==False:
                                Somme=Somme+ViewingAngleZen[ki][kj]
                                cpt=cpt+1
                    if cpt>0:
                        ViewingAngleZen2[i][j]=Somme/(cpt*1.0)



    return [TimeStamp,ViewingAngleZen2,ViewingAngleAzi,float(Zenith),float(Azimuth)]

def RechercheRayleighDat(ThetasS,ThetavS,AziS,AziV):
    Ray_array = [0.23537, 0.16516, 0.09153, 0.04546, 0.03601, 0.0294, 0.02333, 0.01897, 0.01097, 0.00244, 0.0013, 0.00038, 0.01575]

    Sun_azi = float(AziS)
    Sun_zen = float(ThetasS)
    Sat_azi = float(AziV)
    Sat_zen = float(ThetavS)
    a_plus = math.acos(1 * (math.cos(Sun_zen * math.pi / 180) * math.cos(float(Sat_zen) * math.pi / 180) - math.sin(Sun_zen * math.pi / 180) * math.sin(float(Sat_zen) * math.pi / 180) * math.cos(fabs((Sun_azi - float(Sat_azi))) * math.pi / 180))) * 180 / math.pi
    a_minus = math.acos(-1 * (math.cos(Sun_zen * math.pi / 180) * math.cos(float(Sat_zen) * math.pi / 180) - math.sin(Sun_zen * math.pi / 180) * math.sin(float(Sat_zen) * math.pi / 180) * math.cos(fabs((Sun_azi - float(Sat_azi))) * math.pi / 180))) * 180 / math.pi
    if (1 * math.sin(Sun_zen * math.pi / 180) > 1):
        a_trans_sun = math.atan(1 * math.sin(Sun_zen * math.pi / 180)) * 180 / math.pi
    else:
        a_trans_sun = math.asin(1 * math.sin(Sun_zen * math.pi / 180)) * 180 / math.pi
    r_sun_sin_minus = math.sin((Sun_zen - a_trans_sun) * math.pi / 180) * math.sin((Sun_zen - a_trans_sun) * math.pi / 180)
    r_sun_sin_plus = math.sin((Sun_zen + a_trans_sun) * math.pi / 180) * math.sin((Sun_zen + a_trans_sun) * math.pi / 180)
    r_sun_tan_minus = math.tan((Sun_zen - a_trans_sun) * math.pi / 180) * math.tan((Sun_zen - a_trans_sun) * math.pi / 180)
    r_sun_tan_plus = math.tan((Sun_zen + a_trans_sun) * math.pi / 180) * math.tan((Sun_zen + a_trans_sun) * math.pi / 180)
    r_sun = 0.5 * (r_sun_sin_minus / r_sun_sin_plus + r_sun_tan_minus / r_sun_tan_plus)
    # --------------------------------------------------------------------------
    if (float(Sat_zen) == 0):
        r_sat = ((1 - 1) * (1 - 1)) / ((1 + 1) * (1 + 1))
    else:
        a_trans_sat = math.asin(1 * math.sin(float(Sat_zen) * math.pi / 180)) * 180 / math.pi
        r_sat_sin_minus = math.sin((float(Sat_zen) - a_trans_sat) * math.pi / 180) * math.sin((float(Sat_zen) - a_trans_sat) * math.pi / 180)
        r_sat_sin_plus = math.sin((float(Sat_zen) + a_trans_sat) * math.pi / 180) * math.sin((float(Sat_zen) + a_trans_sat) * math.pi / 180)
        r_sat_tan_minus = math.tan((float(Sat_zen) - a_trans_sat) * math.pi / 180) * math.tan((float(Sat_zen) - a_trans_sat) * math.pi / 180)
        r_sat_tan_plus = math.tan((float(Sat_zen) + a_trans_sat) * math.pi / 180) * math.tan((float(Sat_zen) + a_trans_sat) * math.pi / 180)
        r_sat = 0.5 * (r_sat_sin_minus / r_sat_sin_plus + r_sat_tan_minus / r_sat_tan_plus)
    P_plus = 0.75 * (1 + math.cos(a_plus * math.pi / 180) * math.cos(a_plus * math.pi / 180))
    P_minus = 0.75 * (1 + math.cos(a_minus * math.pi / 180) * math.cos(a_minus * math.pi / 180))
    p_ray = P_minus + (r_sun + r_sat) * P_plus
    sub_ray = 4 * math.cos(Sun_zen * math.pi / 180) * math.cos(float(Sat_zen) * math.pi / 180)


    for i in range(len(Ray_array)):
        Ray_array[i] = Ray_array[i] * p_ray / sub_ray
    return Ray_array


#Module de lecture des metadonnees
def LireXML(filename):
    fic=open("%s" %filename,"r")
    lineF = fic.readline().strip()
    cpt=0
    bool=0
    while bool<2:
        lineF = fic.readline().strip()
        if "SENSING_TIME" in lineF:
            bool=bool+1
            TimeStamp=(((lineF.split(">"))[1]).split("<"))[0]
        if "Mean_Viewing_Incidence_Angle bandId=" in lineF:
            bool=bool+1
            lineF = fic.readline().strip()
            Zenith=(((lineF.split(">"))[1]).split("<"))[0]
            lineF = fic.readline().strip()
            Azimuth=(((lineF.split(">"))[1]).split("<"))[0]
        cpt=cpt+1
    fic.close()
    return [TimeStamp,float(Zenith),float(Azimuth)]


def pipe(args_a, args_b):
    p1 = subprocess.Popen(args_a, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(args_b, stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    return p2

def check_output(process):
    """boilerplate for subprocess call"""
    output, error = process.communicate()
    if error is not None and len(error) > 1:
        logging.error("The command `{' '.join(process.args)}` failed:\n", error.decode('utf-8'))
        sys.exit(1)
    if output is not None:
        return output.decode('utf-8')
    logging.error('Error: the command `{" ".join(process.args)}` generated no output')
    sys.exit(1)

###### DEBUT DU PROGRAMME MAIN
if __name__ == "__main__":
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    start_time = time.time()
    logging.info("Progress[%]: 0")
    #Lecture des arguments
    Dir=sys.argv[1]
    TextWrite=sys.argv[2]
    OutPath =sys.argv[4]

    LSB01=glob.glob(os.path.join(Dir, "IMG_DATA", "*B01.jp2"))
    DirW=os.path.dirname(LSB01[0])
    ID=(os.path.basename(LSB01[0])).replace("_B01.jp2","")

    newfolder = OutPath #Dir+"/output_DATA"

    if os.path.isdir(newfolder):
        shutil.rmtree(newfolder)

    os.makedirs(newfolder)
    subprocess.call(['chmod', '-R', '+wrx', newfolder])

    Pressure=getPressure(Dir)

    TextResolution="_60_10_10_10_20_20_20_10_60_60_20_20_20"

    Bande1=None
    Bande2=None
    Bande3=None
    Bande4=None
    Bande5=None
    Bande6=None
    Bande7=None
    Bande8=None
    Bande9=None
    Bande10=None
    Bande11=None
    Bande12=None
    Bande8A=None


    Lamb=[]
    Lamb.append("")
    Lamb.append(0.4439)
    Lamb.append(0.4966)
    Lamb.append(0.56)
    Lamb.append(0.6645)
    Lamb.append(0.7039)
    Lamb.append(0.7402)
    Lamb.append(0.7825)
    Lamb.append(0.8351)
    Lamb.append(0.945)
    Lamb.append(1.3735)
    Lamb.append(1.6137)
    Lamb.append(2.2024)
    Lamb.append(0.8648)

    boolCloud=0

    #Lecture des bandes 1 et 10 pour faire un pseudo masque nuage
    if os.path.isfile(os.path.join(DirW, f"{ID}_B01.jp2")):
        boolCloud=1
        CloudGlob=[]
        CloudGeo=[]
        if os.path.isfile(os.path.join(DirW, f"{ID}_B01.TIF")):
            ds = gdal.Open(os.path.join(DirW, f"{ID}_B01.TIF"), gdal.GA_ReadOnly)

        else:
            ds = gdal.Open(os.path.join(DirW, f"{ID}_B01.jp2"), gdal.GA_ReadOnly)

        band = ds.GetRasterBand(1)
        logging.debug("Temps Avant %s seconds ---" % (time.time() - start_time))
        CloudGlob.append(band.ReadAsArray().astype(np.int))
        CloudGeo.append(ds.GetGeoTransform())



    prj=ds.GetProjection()
    tmp=prj.split("\"")
    EPSG=tmp[len(tmp)-2]
    ds=None



    RepXml=getAngles2(os.path.join(Dir, "MTD_TL.xml"))
    TimeStamp=RepXml[0]
    ArrayZenithV=RepXml[1]
    ArrayAzimuthV=RepXml[2]
    ZenithV=RepXml[3]
    AzimuthV=RepXml[4]
    logging.debug("Temps Apres %s seconds ---" % (time.time() - start_time))




    #Creation des tableaux de Zenith et Azimuth solar angles
    TimeStamp=TimeStamp.split("T")
    Date=TimeStamp[0]
    Heure=(TimeStamp[1].split("."))[0]
    Date=Date.split("-")
    Heure=Heure.split(":")

    TimeStamp=datetime.datetime(int(Date[0]),int(Date[1]),int(Date[2]),int(Heure[0]),int(Heure[1]),int(Heure[2]))

    gdaltrans_command = ['gdaltransform', '-s_srs', f'EPSG:{EPSG}', '-t_srs', 'EPSG:4326']

    ArrayLat=[[0,0],[0,0]]
    ArrayLon=[[0,0],[0,0]]
    Coord=check_output(pipe(['echo', str(CloudGeo[0][0]), str(CloudGeo[0][3])], gdaltrans_command))
    Coord=Coord.split(" ")
    X=Coord[0]
    Y=Coord[1]
    logging.debug("{}\t{}\t{}\t{}".format(X[2:] , Y , type(X) , type(Y)))
    X = X[2:]
    ArrayLon[0][0]=float(X)
    ArrayLat[0][0]=float(Y)

    Coord=check_output(pipe(['echo', str(CloudGeo[0][0]), str(CloudGeo[0][3]+CloudGeo[0][5]*len(CloudGlob[0]))], gdaltrans_command))
    Coord=Coord.split(" ")
    X=Coord[0]
    X = X[2:]
    Y=Coord[1]
    ArrayLon[1][0]=float(X)
    ArrayLat[1][0]=float(Y)

    Coord=check_output(pipe(['echo', str(CloudGeo[0][0]+CloudGeo[0][1]*len(CloudGlob[0][0])), str(CloudGeo[0][3]+CloudGeo[0][5]*len(CloudGlob[0]))], gdaltrans_command))
    Coord=Coord.split(" ")
    X=Coord[0]
    X = X[2:]
    Y=Coord[1]
    ArrayLon[1][1]=float(X)
    ArrayLat[1][1]=float(Y)

    Coord=check_output(pipe(['echo', str(CloudGeo[0][0]+CloudGeo[0][1]*len(CloudGlob[0][0])), str(CloudGeo[0][3])], gdaltrans_command))
    Coord=Coord.split(" ")
    X=Coord[0]
    X = X[2:]
    Y=Coord[1]
    ArrayLon[0][1]=float(X)
    ArrayLat[0][1]=float(Y)

    logging.debug ("{}".format(ArrayLon))
    logging.debug ("{}".format(ArrayLat))
    ArrayLon2=scipy.ndimage.zoom(ArrayLon, 10, order=1)
    ArrayLat2=scipy.ndimage.zoom(ArrayLat, 10, order=1)


    ArrayZenithS=[]
    ArrayRayleigh=[]
    for i in range(0,20):
        tmp=[]
        tmp2=[]
        for j in range(0,20):
            tmp.append([])
            tmp2.append([])
        ArrayZenithS.append(tmp)
        ArrayRayleigh.append(tmp2)

    ArrayAzimuthS=copy.deepcopy(ArrayZenithS)

    costeta=np.zeros((20,20))
    costetav=np.zeros((20,20))

    logging.debug ("{} {}".format(TimeStamp , type(TimeStamp)))
    for i in range(0,20):
        for j in range(0,20):
            ElevationS=get_altitude(ArrayLat2[i][j],ArrayLon2[i][j], TimeStamp.replace(tzinfo=pytz.UTC))
            AzimuthS=get_azimuth(ArrayLat2[i][j],ArrayLon2[i][j], TimeStamp.replace(tzinfo=pytz.UTC))
            ZenithS=90-ElevationS
            ArrayZenithS[i][j]=ZenithS

            AzimuthS=-AzimuthS-180
            ArrayAzimuthS[i][j]=AzimuthS

            DiffAz=(AzimuthS-AzimuthV)*math.pi/180
            deltaphi=math.atan2((1-(math.cos(DiffAz))**2)**0.5,(math.cos(DiffAz)))*180/math.pi
            if np.isnan(ArrayZenithV[i][j])==False:
                #Ray=RechercheRayleigh(ZenithS,ArrayZenithV[i][j],deltaphi)
                Ray=RechercheRayleighDat(ZenithS,ArrayZenithV[i][j],AzimuthS,AzimuthV)
                Rayleigh=['']
                Rayleigh.extend(Ray)
                ArrayRayleigh[i][j]=Rayleigh
                costetav[i][j]=math.cos(ArrayZenithV[i][j]*math.pi/180)
            else:
                ArrayRayleigh[i][j]=['','0','0','0','0','0','0','0','0','0','0','0','0','0']
            costeta[i][j]=math.cos(ZenithS*math.pi/180)


    ArrayRayleigh2=np.array(ArrayRayleigh)

    logging.debug ("{}".format(ArrayRayleigh2[0][0]))
    Resolution=TextResolution.split("_")
    Write=TextWrite.split("_")



    ds = None

    #Lecture des fichiers d'entree indiques pour la correction
    dataGlob=[]
    Rayleigh=[]
    GeoGlob=[]

    logging.debug("Temps Lecture %s seconds ---" % (time.time() - start_time))

    i=1
    while i<=13:
        Resolution[i]=int(Resolution[i])
        Lamb[i]=float(Lamb[i])
        i=i+1


    logging.debug("Temps costeta20 Avant %s seconds ---" % (time.time() - start_time))
    costeta20=scipy.ndimage.zoom(costeta, 274.5, order=3)
    logging.debug("Temps costeta20 Avant %s seconds ---" % (time.time() - start_time))
    costetav20=scipy.ndimage.zoom(costetav, 274.5, order=3)
    logging.debug("Temps costetav20 Avant %s seconds ---" % (time.time() - start_time))

    logging.debug("Temps costeta60 Avant %s seconds ---" % (time.time() - start_time))
    costeta60=scipy.ndimage.zoom(costeta, 91.5, order=3)
    logging.debug("Temps costeta60 Avant %s seconds ---" % (time.time() - start_time))
    costetav60=scipy.ndimage.zoom(costetav, 91.5, order=3)
    logging.debug("Temps costetav60 Avant %s seconds ---" % (time.time() - start_time))


    Moy_O3=max(0,getO3(Dir))
    logging.info ("Ozone moyen : %s DU (avant 332.8 DU) " %Moy_O3)
    Moy_O3=Moy_O3*0.001

    RatioIn=10000
    RatioOut=1.0

    ####################################
    # ECRIRE ICI LE NUM DE LA VERSION DE L ALGO (%BANDES UTILISEES POUR EPSMAX /COR ATM)
    version="_v3"
    ####################################

    #
    costeta10=scipy.ndimage.zoom(costeta20, 2, order=1)
    costetav10=scipy.ndimage.zoom(costetav20, 2, order=1)


    #NoData=np.minimum(np.minimum(np.minimum(dataGlob[Bande2],dataGlob[Bande3]),dataGlob[Bande4]),dataGlob[Bande8])


    #Application des corrections en appliquant le Dark Pixel
    dataRef=np.zeros((3,10980,10980), 'float')
    iBande=0
    # TODO: Note that dataGlob and Rayleigh are not initialized!
    stepNumber = len(Write)
    progress = 0.0
    logging.info("Progress[%]: 5")
    while iBande < stepNumber:
            dataToWrite=[]
            RayToWrite=[]
            CodeBande=Write[iBande]
            if Write[iBande]=="8A":
                    NumBande=int(13)
            else:
                    NumBande=int(Write[iBande])
            if NumBande==1 and ( Bande1 is not None):
                    dataToWrite=copy.deepcopy(dataGlob[Bande1])
                    RayToWrite=copy.deepcopy(Rayleigh[Bande1])
            if NumBande==2 and (Bande2 is not None):
                    dataToWrite=copy.deepcopy(dataGlob[Bande2])
                    RayToWrite=copy.deepcopy(Rayleigh[Bande2])
            if NumBande==3 and ( Bande3 is not None ):
                    dataToWrite=copy.deepcopy(dataGlob[Bande3])
                    RayToWrite=copy.deepcopy(Rayleigh[Bande3])
            if NumBande==4 and ( Bande4 is not None ):
                    dataToWrite=copy.deepcopy(dataGlob[Bande4])
                    RayToWrite=copy.deepcopy(Rayleigh[Bande4])
            if NumBande==5 and ( Bande5 is not None ):
                    dataToWrite=copy.deepcopy(dataGlob[Bande5])
                    RayToWrite=copy.deepcopy(Rayleigh[Bande5])
            if NumBande==6 and ( Bande6 is not None ):
                    dataToWrite=copy.deepcopy(dataGlob[Bande6])
                    RayToWrite=copy.deepcopy(Rayleigh[Bande6])
            if NumBande==7 and ( Bande7 is not None ):
                    dataToWrite=copy.deepcopy(dataGlob[Bande7])
                    RayToWrite=copy.deepcopy(Rayleigh[Bande7])
            if NumBande==8 and ( Bande8 is not None ):
                    dataToWrite=copy.deepcopy(dataGlob[Bande8])
                    RayToWrite=copy.deepcopy(Rayleigh[Bande8])
            if NumBande==9 and ( Bande9 is not None ):
                    dataToWrite=copy.deepcopy(dataGlob[Bande9])
                    RayToWrite=copy.deepcopy(Rayleigh[Bande9])
            if NumBande==10 and ( Bande10 is not None ):
                    dataToWrite=copy.deepcopy(dataGlob[Bande10])
                    RayToWrite=copy.deepcopy(Rayleigh[Bande10])
            if NumBande==11 and ( Bande11 is not None ):
                    dataToWrite=copy.deepcopy(dataGlob[Bande11])
                    RayToWrite=copy.deepcopy(Rayleigh[Bande11])
            if NumBande==12 and ( Bande12 is not None ):
                    dataToWrite=copy.deepcopy(dataGlob[Bande12])
                    RayToWrite=copy.deepcopy(Rayleigh[Bande12])
            if NumBande==13 and ( Bande8A is not None ):
                    dataToWrite=copy.deepcopy(dataGlob[Bande8A])
                    RayToWrite=copy.deepcopy(Rayleigh[Bande8A])
            if dataToWrite == []:
                    if os.path.isfile(os.path.join(DirW, f"{ID}_B{Write[iBande]}.TIF")):
                            ds = gdal.Open(os.path.join(DirW, f"{ID}_B{Write[iBande]}.TIF"), gdal.GA_ReadOnly)
                    else:
                            ds = gdal.Open(os.path.join(DirW, f"{ID}_B{Write[iBande]}.jp2"), gdal.GA_ReadOnly)

                    logging.info(f"Reading Band {Write[iBande]}")
                    band = ds.GetRasterBand(1)
                    SizeX = ds.RasterXSize
                    logging.debug("Temps Avant %s seconds ---" % (time.time() - start_time))
                    dataToWrite=band.ReadAsArray()
                    logging.debug("Temps Apres %s seconds ---" % (time.time() - start_time))
                    ExtRay=ArrayRayleigh2[...,int(NumBande)].astype(np.float)
                    RayToWrite=scipy.ndimage.zoom(ExtRay, SizeX/20.0, order=3)
            logging.debug("Temps Start %s seconds ---" % (time.time() - start_time))
            tauR=8.524*10**(-3)*Lamb[NumBande]**(-4)+9.63*10**(-5)*Lamb[NumBande]**(-6)+1.1*10**(-6)*Lamb[NumBande]**(-8)
            tau_O3=getTau(Lamb[NumBande])
            logging.info ("Ray %s : %s" %(NumBande,RayToWrite[0][0]))
            if int(Resolution[NumBande])==10:
                    Correc=(dataToWrite/(RatioIn*np.exp(-Moy_O3*tau_O3*(1/costeta10+1/costetav10)))-RayToWrite)/(np.exp(-0.5*tauR*(1/costeta10+1/costetav10)))*RatioOut
                    logging.debug("Temps Calc %s seconds ---" % (time.time() - start_time))
                    #Correc[NoData<=0]=-99
            elif int(Resolution[NumBande])==20:
                    Correc=(dataToWrite/(RatioIn*np.exp(-Moy_O3*tau_O3*(1/costeta20+1/costetav20)))-RayToWrite)/(np.exp(-0.5*tauR*(1/costeta20+1/costetav20)))*RatioOut
                    #Correc[dataToWrite<=0]=-99
            elif int(Resolution[NumBande])==60:
                    Correc=(dataToWrite/(RatioIn*np.exp(-Moy_O3*tau_O3*(1/costeta60+1/costetav60)))-RayToWrite)/(np.exp(-0.5*tauR*(1/costeta60+1/costetav60)))*RatioOut
                    #Correc[dataToWrite<=0]=-99
            #Enregistrement des fichiers corriges - surface reflectances
            driver = gdal.GetDriverByName('GTiff')
            inDs = gdal.Open(os.path.join(DirW, f"{ID}_B{Write[iBande]}.jp2"))
            outDs = driver.Create(os.path.join(newfolder, f"C_{ID}_B{Write[iBande]}.TIF"), len(Correc[1]), len(Correc), 1, gdal.GDT_Float32)
            outBand = outDs.GetRasterBand(1)
            outBand.WriteArray(Correc)
            outBand.FlushCache()
            outBand.SetNoDataValue(-99)
            outDs.SetGeoTransform(inDs.GetGeoTransform())
            outDs.SetProjection(inDs.GetProjection())
            logging.info ("Band %s SAVED" %Write[iBande])
            #Correc = Correc.astype(int)
            if int(NumBande)==2:
                    dataRef[0] = Correc
            if int(NumBande)==3:
                    dataRef[1] = Correc
            if int(NumBande)==4:
                    dataRef[2] = Correc
            currentprogress = 90*iBande/float(stepNumber)
            if(currentprogress>progress):
                progress=int(currentprogress)
                current_progress = (progress+5)
                logging.info("Progress[%]: {}".format(current_progress))
            iBande+=1
    #Creation du jpg de l'image corrigee (BOA)
    PixBlanc=float(sys.argv[3])

    rgbArray = np.zeros((len(dataRef[0]),len(dataRef[0][1]),3), 'uint8')
    rgbArray[..., 0] = np.maximum(0,np.minimum(dataRef[2]*255/(RatioOut*PixBlanc),255))  #Rouge
    rgbArray[..., 1] = np.maximum(0,np.minimum(dataRef[1]*255/(RatioOut*PixBlanc),255))  #Vert
    rgbArray[..., 2] = np.maximum(0,np.minimum(dataRef[0]*255/(RatioOut*PixBlanc),255))  #Bleu
    img = Image.fromarray(rgbArray)
    img.save(os.path.join(newfolder, f'RGB_{ID}.jpg'), quality=80)

    RGB = np.zeros((3,len(dataRef[0]),len(dataRef[0][1])), 'uint16')
    RGB[0] = np.maximum(0,np.minimum(dataRef[2]*255/(RatioOut*PixBlanc),255))  #Rouge
    RGB[1] = np.maximum(0,np.minimum(dataRef[1]*255/(RatioOut*PixBlanc),255))  #Vert
    RGB[2] = np.maximum(0,np.minimum(dataRef[0]*255/(RatioOut*PixBlanc),255))  #Bleu


    inDs = gdal.Open(os.path.join(DirW, f"{ID}_B02.jp2"))
    inGs=inDs.GetGeoTransform()

    outDs = driver.Create(os.path.join(newfolder, f"RGB_{ID}.TIF"), len(dataRef[0][1]), len(dataRef[0]), 3, gdal.GDT_Byte)
    outDs.GetRasterBand(1).WriteArray(RGB[0])
    outDs.GetRasterBand(2).WriteArray(RGB[1])
    outDs.GetRasterBand(3).WriteArray(RGB[2])
    outDs.SetGeoTransform(inGs)
    outDs.SetProjection(inDs.GetProjection())
    outDs.FlushCache()



    try:
        DimTile=109800.0/(len(ArrayRayleigh2)*1.0)
        #Enregistrement du Eps et MinB11
        inDs = gdal.Open(os.path.join(DirW, f"{ID}_B01.jp2"))
        inGs=inDs.GetGeoTransform()

        GeoT=[inGs[0],DimTile,inGs[2],inGs[3],inGs[4],-DimTile]

        outDs = driver.Create(os.path.join(newfolder, f"Rayleigh865_{ID}.TIF"), len(ArrayRayleigh2[1]), len(ArrayRayleigh2), 1, gdal.GDT_Float32, options=['TILED=YES','BIGTIFF=IF_SAFER','BLOCKXSIZE=512','BLOCKYSIZE=512','COMPRESS=LZW'])
        outBand = outDs.GetRasterBand(1)
        outBand.WriteArray(ArrayRayleigh2[..., 13])
        outBand.FlushCache()
        outBand.SetNoDataValue(-99)
        outDs.SetGeoTransform(GeoT)
        outDs.SetProjection(inDs.GetProjection())

        logging.info ("Rayleigh SAVED")
        logging.info("Progress[%]: 96")


        DimTile=109800.0/(len(ArrayZenithS)*1.0)
        GeoT=[inGs[0],DimTile,inGs[2],inGs[3],inGs[4],-DimTile]


        driver.Create(os.path.join(newfolder, f"SolarZenith_{ID}.TIF"), len(ArrayZenithS[1]), len(ArrayZenithS), 1, gdal.GDT_Float32, options=['TILED=YES','BIGTIFF=IF_SAFER','BLOCKXSIZE=512','BLOCKYSIZE=512','COMPRESS=LZW'])
        outBand = outDs.GetRasterBand(1)
        outBand.WriteArray(np.array(ArrayZenithS))
        outBand.FlushCache()
        outBand.SetNoDataValue(-99)
        outDs.SetGeoTransform(GeoT)
        outDs.SetProjection(inDs.GetProjection())
        logging.info ("ArrayZenithS SAVED")
        logging.info("Progress[%]: 97")
        DimTile=109800.0/(len(ArrayAzimuthS)*1.0)
        GeoT=[inGs[0],DimTile,inGs[2],inGs[3],inGs[4],-DimTile]
        outDs = driver.Create(os.path.join(newfolder, f"SolarAzimuth_{ID}.TIF"), len(ArrayAzimuthS[1]), len(ArrayAzimuthS), 1, gdal.GDT_Float32, options=['TILED=YES','BIGTIFF=IF_SAFER','BLOCKXSIZE=512','BLOCKYSIZE=512','COMPRESS=LZW'])
        outBand = outDs.GetRasterBand(1)
        outBand.WriteArray(np.array(ArrayAzimuthS))
        outBand.FlushCache()
        outBand.SetNoDataValue(-99)
        outDs.SetGeoTransform(GeoT)
        outDs.SetProjection(inDs.GetProjection())
        logging.info ("SolarZenith SAVED")
        logging.info("Progress[%]: 98")
        DimTile=109800.0/(len(ArrayZenithV)*1.0)
        GeoT=[inGs[0],DimTile,inGs[2],inGs[3],inGs[4],-DimTile]
        outDs = driver.Create(os.path.join(newfolder, f"ViewingZenith_{ID}.TIF"), len(ArrayZenithV[1]), len(ArrayZenithV), 1, gdal.GDT_Float32, options=['TILED=YES','BIGTIFF=IF_SAFER','BLOCKXSIZE=512','BLOCKYSIZE=512','COMPRESS=LZW'])
        outBand = outDs.GetRasterBand(1)
        outBand.WriteArray(np.array(ArrayZenithV))
        outBand.FlushCache()
        outBand.SetNoDataValue(-99)
        outDs.SetGeoTransform(GeoT)
        outDs.SetProjection(inDs.GetProjection())
        logging.info ("SolarZenith SAVED")
        logging.info("Progress[%]: 99")
        DimTile=109800.0/(1.0*1.0)
        GeoT=[inGs[0],DimTile,inGs[2],inGs[3],inGs[4],-DimTile]
        outDs = driver.Create(os.path.join(newfolder, f"ViewingAzimuth_{ID}.TIF"), 1, 1, 1, gdal.GDT_Float32, options=['TILED=YES','BIGTIFF=IF_SAFER','BLOCKXSIZE=512','BLOCKYSIZE=512','COMPRESS=LZW'])
        outBand = outDs.GetRasterBand(1)
        outBand.WriteArray(np.array([[AzimuthV]]))
        outBand.FlushCache()
        outBand.SetNoDataValue(-99)
        outDs.SetGeoTransform(GeoT)
        outDs.SetProjection(inDs.GetProjection())
        logging.info ("ViewingAzimuth SAVED")
    except:
        logging.info ("Warning: not enough information")
    logging.info("Progress[%]: 100")
    sys.exit(0)