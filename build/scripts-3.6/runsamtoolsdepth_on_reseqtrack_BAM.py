#!/nfs/production/reseq-info/work/ernesto/bin/anaconda3/envs/ehive/bin/python
'''
Created on 12 Oct 2016

This script is used to run samtools depth on a file that is in a Reseqtrack DB

@author: ernesto
'''
from BamQC import BamQC,SDepth
from pprint import pprint
from ReseqTrackDB import *
import argparse

#get command line arguments

#RESEQTRACK DB conn params
parser = argparse.ArgumentParser(description='Script to calculate different QC metrics on the WGS BAMs')
parser.add_argument('--hostname', type=str, required=True, help='Hostname for ReseqTrack DB' )
parser.add_argument('--username', type=str, required=True, help='User for ReseqTrack DB' )
parser.add_argument('--port', type=int, required=True, help='Port number in the ReseqTrack DB' )
parser.add_argument('--pwd', type=str, help='PWD for the ReseqTrack DB' )
parser.add_argument('--db', type=str, required=True, help='DB name in the ReseqTrack DB' )

parser.add_argument('--samtools', type=str, required=True, help='Folder containing the SAMtools binary' )
parser.add_argument('--filename', type=str, required=True, help='Base name of file that will be analysed' )
parser.add_argument('--list', type=argparse.FileType('r'), required=True, help='File with a list of chros to be analysed (comma-separated)' )


args = parser.parse_args()

if __name__ == '__main__':
    
    hostname=args.hostname
    username=args.username
    db=args.db
    port=args.port
    pwd=args.pwd
    
    reseqdb = ReseqTrackDB(host=hostname,user=username,port=port,pwd=pwd,db=db)
    
    file=reseqdb.fetch_file_by_filename(args.filename)
    
    #list of chros to analyse
    list_chros=[]
    for line in args.list:
        line=line.strip("\n")
        list_chros=line.split(',')
    
    bam = BamQC(file.path,args.samtools)
    
    list_of_samples=bam.list_of_samples()
    list_of_readgroups=bam.list_of_readgroups()
    stats=bam.get_simple_stats()
    for k,v in stats.items(): Attribute(table_name="file",other_id=file.dbID,name=k,value=v).store(reseqdb)
    for i in range(len(list_of_readgroups)): Attribute(table_name="file",other_id=file.dbID,name="readgroup%d" % (i+1),value=list_of_readgroups[i]).store(reseqdb)
    for i in range(len(list_of_samples)): Attribute(table_name="file",other_id=file.dbID,name="sample%d" % (i+1),value=list_of_samples[i]).store(reseqdb)
    cov_list=bam.run_samtools_depth(chros=list_chros)
    aggr_cov=bam.aggregate_stats(cov_list)
    Attribute(table_name="file",other_id=file.dbID,name="COV:contig",value=aggr_cov.contig).store(reseqdb)
    Attribute(table_name="file",other_id=file.dbID,name="COV:breadth",value=aggr_cov.breadth).store(reseqdb)
    Attribute(table_name="file",other_id=file.dbID,name="COV:depth",value=aggr_cov.depth).store(reseqdb)
    Attribute(table_name="file",other_id=file.dbID,name="COV:max",value=aggr_cov.max).store(reseqdb)
    
    pass
