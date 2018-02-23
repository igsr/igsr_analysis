'''
Utility to create a new entry in the Reseqtrackdb's Collection table that will be associated to a certain type of files
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

parser.add_argument('--ftype_1', type=str, required=True, help='File type1 that will be used to generate the new collection' )
parser.add_argument('--ftype_2', type=str, required=True, help='File type2 that will be used to generate the new collection' )
parser.add_argument('--ftype_3', type=str, required=False, help='File type3 that will be used to generate the new collection' )
parser.add_argument('--collection_type', type=str, required=True, help='Type for the new collection' )
parser.add_argument('--collection_name', type=str, required=True, help='Name for the new collection' )


args = parser.parse_args()

def get_files_by_type(reseqdb,type):
    '''
    
    Parameters
    ----------
    reseqdb : ReseqTrackDB object, Required
    type : str, Required
           File type to retrieve

    Returns
    -------
    A list with dbIDs of the files

    '''

    l=reseqdb.fetch_files_by_type(type)

    return [x.dbID for x in l]

if __name__ == '__main__':
    log_filename="create_collection_to_transpose.log"

    logger = logging.getLogger("col_2_transpose")
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
    
    l1=get_files_by_type(reseqdb,args.ftype_1)
    l2=get_files_by_type(reseqdb,args.ftype_2)
    
    l3=[]
    if args.ftype_3:
        l3=get_files_by_type(reseqdb,args.ftype_3)
    
    others_ids=list(set(l1+l2+l3))

    new_c=Collection(name=args.collection_name,type=args.collection_type,others_dbIDs=others_ids,table_name='file')

    new_c.store(reseqdb)

    logger.info("Done")
