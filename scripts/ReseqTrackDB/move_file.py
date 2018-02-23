'''
Utility to move a file to a new location, and this modified file will be stored in the database
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

parser.add_argument('--filepath', type=str, required=True, help='Filepath for file to be moved' )
parser.add_argument('--newpath', type=str, required=True, help='New path for file' )
parser.add_argument('--type', type=str, required=True, help='Type for file' )


args = parser.parse_args()

if __name__ == '__main__':
    log_filename="move_file.log" 

    logger = logging.getLogger("movefile")
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
    f=File(path=args.filepath,type=args.type)

    f.move(reseqdb, newpath= args.newpath, do_md5=True)

    logger.info("Done")
