#python LAC_S2.py -h to see description 

import os, sys, argparse,subprocess,glob
from datetime import datetime,timedelta


parser = argparse.ArgumentParser(description='Download and process S2 data')
parser.add_argument('-dirout','--dirout', required=True, help='(e.g: /home)')
parser.add_argument('-atmcor','--atmcor', required=True, help='(e.g: y means full correction while n means only Rayleigh)')
parser.add_argument('-bands','--bands', required=False, help='(e.g: Set the bands you want as output : 02_03_04_11 as default)')

args = parser.parse_args()

DirOut=str(args.dirout)
Bands = str(args.bands)

Atmcor=str(args.atmcor)
# Level=str(args.level)

Level="L1C"
PixSat="0.15"
#DirScript="/mount/internal/work-st/projects/esrin-079/1402-seom/S2/scripts/S2_L8/S2"
DirScript=os.path.dirname(os.path.realpath(__file__))
		
#print "python "+DirScript+"/Sentinel_download2.py --lat "+Lat+" --lon "+Lon+" -a "+DirScript+"/apihub.txt -t "+Code+" -l "+Level+" -d "+TDateTemp+" -f "+TDateTemp+" -w "+DirOut+"/"+Code+"/"+TDateTemp+" --dhus"
#test=subprocess.check_output("python "+DirScript+"/Sentinel_download2.py --lat "+Lat+" --lon "+Lon+" -a "+DirScript+"/apihub.txt -t "+Code+" -l "+Level+" -d "+TDateTemp+" -f "+TDateTemp+" -w "+DirOut+"/"+Code+"/"+TDateTemp+" --dhus", shell=True)
		
if Atmcor=="n":
	LS=glob.glob(DirOut+"/MTD_TL.xml")
	if len(LS)>0: #TEST pour voir s'il y a un fichier telecharge a la date donnee
		#test=subprocess.check_output("python "+DirScript+"/SEOM_Rayleigh_SCIHUB.py "+DirOut+"/"+Code+"/"+TDateTemp+" \""+Bands+"\" "+PixSat, shell=True)
		print ("python "+DirScript+"/SEOM_Rayleigh_LAC.py "+DirOut+" \""+Bands+"\" "+PixSat)
		test=subprocess.Popen("python "+DirScript+"/SEOM_Rayleigh_LAC.py "+DirOut+" \""+Bands+"\" "+PixSat,stderr=subprocess.PIPE,shell=True)
		output,error=test.communicate()
		print output
else:
	LS=glob.glob(DirOut+"/MTD_TL.xml")
	if len(LS)>0: #TEST pour voir s'il y a un fichier telecharge a la date donnee
		#test=subprocess.check_output("python "+DirScript+"/SEOM_Rayleigh_SCIHUB.py "+DirOut+"/"+Code+"/"+TDateTemp+" \""+Bands+"\" "+PixSat, shell=True)
		test=subprocess.Popen("python "+DirScript+"/SEOM_RayleighAerosols_LAC.py "+DirOut+" \""+Bands+"\" 02_03_04_08",stderr=subprocess.PIPE,shell=True)
		output,error=test.communicate()
		print output
