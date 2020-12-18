#python LAC_S2.py -h to see description 

import os, sys, argparse, subprocess
from datetime import datetime,timedelta
import presteps, poststeps

def parse_args(argv, struct):
    args = {}
    argv = [arg.replace('-', '').replace(' ', '').strip() if arg.startswith('-') else arg for arg in argv]

    for key in struct:
        if key in argv:
            i = argv.index(key)
            if i+1 < len(argv) and argv[i+1] not in struct.keys():
                args[key] = argv[i+1]
            else:
                print(f'Misformed argument: -{key} or --{key}')
                sys.exit(1)
        else:
            if struct[key][1] is None:
                print(f'Missing argument: -{key} or --{key}')
                sys.exit(1)
            else:
                args[key] = struct[key][1]
    return args


argv = sys.argv[1:]

arg_struct = {
    'dirin': ('(e.g: /home)', None),
    'dirout': ('e.g: /home)', None),
    'atmcor': ('(e.g: y means full correction while n means only Rayleigh)', None),
    'bands': ('(e.g: Set the bands you want as output : 02_03_04_11 as default)', '02_03_04_11'),
}



args = parse_args(argv, arg_struct)

DirIn =str(args['dirin'])
Bands = str(args['bands'])
Atmcor=str(args['atmcor'])
DirOut=str(args['dirout'])

presteps.prepare(DirIn, DirOut, Bands, Atmcor)

TmpDirOut = os.path.join(DirOut, 'TMP')

# Level=str(args.level)

print('Ready to process')
Level="L1C"
PixSat="0.15"
#DirScriptDirScript="/mount/internal/work-st/projects/esrin-079/1402-seom/S2/scripts/S2_L8/S2"
		
#print "python "+DirScript+"/Sentinel_download2.py --lat "+Lat+" --lon "+Lon+" -a "+DirScript+"/apihub.txt -t "+Code+" -l "+Level+" -d "+TDateTemp+" -f "+TDateTemp+" -w "+DirOut+"/"+Code+"/"+TDateTemp+" --dhus"
#test=subprocess.check_output("python "+DirScript+"/Sentinel_download2.py --lat "+Lat+" --lon "+Lon+" -a "+DirScript+"/apihub.txt -t "+Code+" -l "+Level+" -d "+TDateTemp+" -f "+TDateTemp+" -w "+DirOut+"/"+Code+"/"+TDateTemp+" --dhus", shell=True)
#print ("python "+" SEOM_Rayleigh_LAC.py "+DirOut+" \""+Bands+"\" "+PixSat)
if os.path.exists(os.path.join(DirIn, "MTD_TL.xml")): #TEST pour voir s'il y a un fichier telecharge a la date donnee
    if Atmcor=="n":	
        print('Process without full correction')
        #test=subprocess.check_output("python "+DirScript+"/SEOM_Rayleigh_SCIHUB.py "+DirOut+"/"+Code+"/"+TDateTemp+" \""+Bands+"\" "+PixSat, shell=True)
        test=subprocess.Popen(["python", "SEOM_Rayleigh_LAC.py", DirIn, Bands, PixSat, TmpDirOut], stderr=subprocess.PIPE)
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
        test=subprocess.Popen(["python", "SEOM_RayleighAerosols_LAC.py", DirIn, Bands, "02_03_04_08_11", TmpDirOut], stderr=subprocess.PIPE) #+"\" 02_03_04_08"
        #print ("python "+" SEOM_RayleighAerosols_LAC.py "+DirIn+" \""+Bands+"\" 02_03_04_08_11 "+DirOut)
        output,error=test.communicate()
        if error is not None:
            error = error.decode('utf-8')
            if 'Warning' not in error:
                print("The process failed:\n", error)
                sys.exit(1)
        if output is not None:
            print (output)
    poststeps.poststeps(DirOut)
else:
    print(f'Error: input directory `{DirIn}` does not contain the `MTD_TL.xml` file')
    sys.exit(1)