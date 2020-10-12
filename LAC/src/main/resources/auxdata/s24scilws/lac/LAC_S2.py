#python LAC_S2.py -h to see description 

import os, sys, argparse, subprocess
from datetime import datetime,timedelta

parser = argparse.ArgumentParser(description='Download and process S2 data')
parser.add_argument('-dirin','--dirin', required=True, help='(e.g: /home)')
parser.add_argument('-atmcor','--atmcor', required=True, help='(e.g: y means full correction while n means only Rayleigh)')
parser.add_argument('-bands','--bands', required=False, help='(e.g: Set the bands you want as output : 02_03_04_11 as default)')
parser.add_argument('-dirout','--dirout', required=True, help='(e.g: /home)')

print("***HERE***")
print(sys.args)

args = parser.parse_args()

DirIn =str(args.dirin)
Bands = str(args.bands)
Atmcor=str(args.atmcor)
DirOut=str(args.dirout)


# Level=str(args.level)

Level="L1C"
PixSat="0.15"
#DirScriptDirScript="/mount/internal/work-st/projects/esrin-079/1402-seom/S2/scripts/S2_L8/S2"
		
#print "python "+DirScript+"/Sentinel_download2.py --lat "+Lat+" --lon "+Lon+" -a "+DirScript+"/apihub.txt -t "+Code+" -l "+Level+" -d "+TDateTemp+" -f "+TDateTemp+" -w "+DirOut+"/"+Code+"/"+TDateTemp+" --dhus"
#test=subprocess.check_output("python "+DirScript+"/Sentinel_download2.py --lat "+Lat+" --lon "+Lon+" -a "+DirScript+"/apihub.txt -t "+Code+" -l "+Level+" -d "+TDateTemp+" -f "+TDateTemp+" -w "+DirOut+"/"+Code+"/"+TDateTemp+" --dhus", shell=True)
#print ("python "+" SEOM_Rayleigh_LAC.py "+DirOut+" \""+Bands+"\" "+PixSat)
if os.path.exists(os.path.join(DirIn, "MTD_TL.xml")): #TEST pour voir s'il y a un fichier telecharge a la date donnee
    if Atmcor=="n":	
        #test=subprocess.check_output("python "+DirScript+"/SEOM_Rayleigh_SCIHUB.py "+DirOut+"/"+Code+"/"+TDateTemp+" \""+Bands+"\" "+PixSat, shell=True)
        test=subprocess.Popen(["python", "SEOM_Rayleigh_LAC.py", DirIn, Bands, PixSat, DirOut], stderr=subprocess.PIPE)
        output,error=test.communicate()
        if error is not None:
            error = error.decode('utf-8')
            if 'Warning' not in error:
                print("The process failed:\n", error)
                sys.exit(1)
        if output is not None:
            print (output)
    elif Atmcor=="y":
        #test=subprocess.check_output("python "+DirScript+"/SEOM_Rayleigh_SCIHUB.py "+DirOut+"/"+Code+"/"+TDateTemp+" \""+Bands+"\" "+PixSat, shell=True)
        test=subprocess.Popen(["python", "SEOM_RayleighAerosols_LAC.py", DirIn, Bands, "02_03_04_08_11", DirOut], stderr=subprocess.PIPE) #+"\" 02_03_04_08"
        #print ("python "+" SEOM_RayleighAerosols_LAC.py "+DirIn+" \""+Bands+"\" 02_03_04_08_11 "+DirOut)
        output,error=test.communicate()
        if error is not None:
            error = error.decode('utf-8')
            if 'Warning' not in error:
                print("The process failed:\n", error)
                sys.exit(1)
        if output is not None:
            print (output)
