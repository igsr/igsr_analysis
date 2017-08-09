'''

This script will fetch from a ReseqTrackDB all sample-level BAM files from a certain type and will compare their 'readgroup' attribute with the RUN_IDs for  
that particular sample extracted from an IGSR Sequence Index. It will complain if there are RUN_IDs present in the Sequence Index that are not in the DB

Created on 30 Nov 2016

@author: ernesto

'''

import argparse
from SequenceIndex import SequenceIndex
from ReseqTrackDB import ReseqTrackDB

#get command line arguments

#RESEQTRACK DB conn params
parser = argparse.ArgumentParser(description='This script will fetch from a ReseqTrackDB all sample-level BAM files from a certain type and will\
compare their readgroup attribute with the RUN_IDs for that particular sample extracted from an IGSR Sequence Index. It will complain if there are\
 RUN_IDs present in the Sequence Index that are not in the DB')
parser.add_argument('--hostname', type=str, required=True, help='Hostname for ReseqTrack DB' )
parser.add_argument('--username', type=str, required=True, help='User for ReseqTrack DB' )
parser.add_argument('--port', type=int, required=True, help='Port number in the ReseqTrack DB' )
parser.add_argument('--pwd', type=str, help='PWD for the ReseqTrack DB' )
parser.add_argument('--db', type=str, required=True, help='DB name in the ReseqTrack DB' )

parser.add_argument('--type', type=str, required=True, help='Type in the file table that will be analysed' )
parser.add_argument('--index', type=str, required=True, help='Sequence index used for the comparison' )
parser.add_argument('--analysis_group', type=str, required=True, help='Subgroup by analysis group in the sequence index' )

args = parser.parse_args()

reseqdb = ReseqTrackDB(host=args.hostname,user=args.username,port=args.port,pwd=args.pwd,db=args.db)

files=reseqdb.fetch_files_by_type(args.type)

s=SequenceIndex(args.index)

data=s.runs_per_sample(args.analysis_group)

for f in files:
    #get sample id, assuming that the sample name is the first bit of the file name
    sample_name=f.name.split('.')[0]
    attbs=reseqdb.fetch_attributes_by_other_id(f.dbID)
    sel_attbs=filter(lambda x:'readgroup' in x.name, attbs)
    #get set with run_ids in DB
    runs_db=set([a.value for a in sel_attbs])
    
    #get set with run_ids for this sample in the index
    runs_index=data.get(sample_name)
    
    if (len(runs_db ^ runs_index)>0):
        print "DB: %s" %runs_db
        print "INDEX: %s" %runs_index
        raise Exception("run_ids discrepancies for %s" %f.name)
    else:      
        print "[INFO] %s - OK"% f.name

