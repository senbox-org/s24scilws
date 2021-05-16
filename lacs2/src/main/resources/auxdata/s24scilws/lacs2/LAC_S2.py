#python LAC_S2.py -h to see description 
import os, sys, argparse, subprocess
from datetime import datetime,timedelta
import presteps, poststeps
import logging

logging.basicConfig()
logger = logging.getLogger('logger')

def parse_args(argv, struct):
    logging.info('Start LAC_S2 cmd')
    args = {}
    argv = [arg.replace('-', '').replace(' ', '').strip() if arg.startswith('-') else arg for arg in argv]

    for key in struct:
        if key in argv:
            i = argv.index(key)
            if i+1 < len(argv) and argv[i+1] not in struct.keys():
                args[key] = argv[i+1]
            else:
                logging.error('Misformed argument: -{key} or --{key}')
                sys.exit(1)
        else:
            if struct[key][1] is None:
                logging.error('Missing argument: -{key} or --{key}')
                sys.exit(1)
            else:
                args[key] = struct[key][1]
    return args

if __name__ == "__main__":
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    python_exe = sys.executable
    argv = sys.argv[1:]

    arg_struct = {
        'dirin': ('(e.g: /home)', None),
        'dirout': ('e.g: /home)', None),
        'atmcor': ('(e.g: y means full correction while n means only Rayleigh)', None),
        'bands': ('(e.g: Set the bands you want as output : 02_03_04_11 as default)', '02_03_04_11'),
    }



    args = parse_args(argv, arg_struct)

    DirIn =str(args['dirin'])
    ProdName = os.path.split(DirIn)[-1]
    Bands = str(args['bands'])
    Atmcor=str(args['atmcor'])
    DirOut= str(args['dirout']) #, f"{ProdName}_{Bands}_{Atmcor}_LAC")

    presteps.prepare(DirIn, DirOut, Bands, Atmcor)

    GrnPath = os.path.join(DirIn, 'GRANULE')
    GrnDirIn = os.path.join(GrnPath, os.listdir(GrnPath)[0])

    TmpDirOut = os.path.join(DirOut, 'TMP')

    # Level=str(args.level)

    logging.info('Ready to process')
    Level="L1C"
    PixSat="0.15"
    if os.path.exists(os.path.join(GrnDirIn, "MTD_TL.xml")):
        if Atmcor=="n":
            logging.info('Process without full correction')
            src_file = sys.argv[0].replace("LAC_S2.py","SEOM_Rayleigh_LAC.py")
            logging.info(src_file)
            cmd = [python_exe, src_file, GrnDirIn, Bands, PixSat, TmpDirOut]
            try:
                subprocess.check_call(cmd,0)
            except Exception as ex:
                logging.error("The process failed:\n {}".format(ex))
                sys.exit(1)
        elif Atmcor=="y":
            src_file = sys.argv[0].replace("LAC_S2.py","SEOM_RayleighAerosols_LAC.py")
            logging.info(src_file)
            cmd = [python_exe, src_file, GrnDirIn, Bands, "02_03_04_08_11", TmpDirOut]
            try:
                subprocess.check_call(cmd,0)
            except Exception as ex:
                logging.error("The process failed:\n {}".format(ex))
                sys.exit(1)
        poststeps.poststeps(DirOut)
    else:
        logging.error('Error: input directory `{DirIn}` does not contain the `MTD_TL.xml` file')
        sys.exit(1)
    logging.info("LAC processing: SUCCESS")