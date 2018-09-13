'''
Utility to change file type in the DB for some files that are passed through a file containing the basenames that want to be modified
'''
from ReseqTrackDB import *
import argparse
import logging
import os

#RESEQTRACK DB conn params
parser = argparse.ArgumentParser(description='Script to calculate different QC metrics on BAMs')
parser.add_argument('--hostname', type=str, required=True, help='Hostname for ReseqTrack DB' )
parser.add_argument('--username', type=str, required=True, help='User for ReseqTrack DB' )
parser.add_argument('--port', type=int, required=True, help='Port number in the ReseqTrack DB' )
parser.add_argument('--pwd', type=str, help='PWD for the ReseqTrack DB' )
parser.add_argument('--db', type=str, required=True, help='DB name in the ReseqTrack DB' )

parser.add_argument('--ifile', type=str, required=True, help='File with filenames to be modified' )
parser.add_argument('--newtype', type=str, required=True, help='New file type' )

args = parser.parse_args()

if __name__ == '__main__':
    log_filename="change_type.log" 

    logger = logging.getLogger("changetype")
    logger.setLevel(logging.INFO)

    #  create the logging file handler
    fh = logging.FileHandler(log_filename)

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)

    # add handler to logger object
    logger.addHandler(fh)

    logger.info("Program started")

    hostname=args.hostname
    username=args.username
    db=args.db
    port=args.port
    pwd=args.pwd

    reseqdb = ReseqTrackDB(host=hostname,user=username,port=port,pwd=pwd,db=db)


    with open(args.ifile) as f:
        for p in f:
            p=p.rstrip("\n")
            filename=os.path.basename(p)
            logger.info("Changing file type for %s" % filename)
            fileO=reseqdb.fetch_file_by_filename(filename)

            old_type=fileO.type
            newtype=""
            
            if old_type=="LC_BAM" or old_type=="EX_BAM":
                newtype=args.newtype
            else:
                newtype=old_type+"_"+args.newtype
            if not newtype: raise Exception("New type is not set. So file type could not be updated")
            fileO.change_file_type(reseqdb,newtype)
                
                

    logger.info("Done")
