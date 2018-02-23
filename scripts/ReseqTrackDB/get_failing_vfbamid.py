'''

This script will get the paths of the files failing the verifybamid QA criteris

Created on 01 Feb 2017

@author: ernesto

'''

import argparse
from ReseqTrackDB import ReseqTrackDB

#get command line arguments

#RESEQTRACK DB conn params
parser = argparse.ArgumentParser(description='This script will create a Boxplot with the median coverage calculated over a certain BAM file')
parser.add_argument('--hostname', type=str, required=True, help='Hostname for ReseqTrack DB' )
parser.add_argument('--username', type=str, required=True, help='User for ReseqTrack DB' )
parser.add_argument('--port', type=int, required=True, help='Port number in the ReseqTrack DB' )
parser.add_argument('--pwd', type=str, help='PWD for the ReseqTrack DB' )
parser.add_argument('--db', type=str, required=True, help='DB name in the ReseqTrack DB' )

parser.add_argument('--type', type=str, required=True, help='Type in the file table that will be analysed' )
parser.add_argument('--chipmix_thr', type=str, required=True, help='chipmix threshold value. E.g. 0.02' )
parser.add_argument('--freemix_thr', type=str, required=True, help='freemix threshold value. E.g. 0.03' )

args = parser.parse_args()

reseqdb = ReseqTrackDB(host=args.hostname,user=args.username,port=args.port,pwd=args.pwd,db=args.db)

files=reseqdb.fetch_files_by_type(args.type)

medians=[]
means=[]

for f in files:
    attbs=reseqdb.fetch_attributes_by_other_id(f.dbID)
    lchip=list(filter(lambda x:'ALL:chipmix' in x.name, attbs))
    lfree=list(filter(lambda x:'ALL:freemix' in x.name, attbs))
    for chip, free in zip(lchip, lfree):
        if chip.value>float(args.chipmix_thr) and free.value>float(args.freemix_thr):
            print(f.name+"\t"+chip.name+"\t"+str(chip.value)+"\t"+free.name+"\t"+str(free.value))
