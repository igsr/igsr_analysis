#!/nfs/production/reseq-info/work/ernesto/bin/anaconda3/envs/ehive/bin/python
'''
Created on 12 Oct 2016

This script is used to run samtools depth on a BAM file

@author: ernesto
'''
from BamQC import BamQC,SDepth
from pprint import pprint

import argparse

#get command line arguments

#RESEQTRACK DB conn params
parser = argparse.ArgumentParser(description='Script to calculate different QC metrics on the WGS BAMs')

parser.add_argument('--samtools', type=str, required=True, help='Folder containing the SAMtools binary' )
parser.add_argument('--filename', type=str, required=True, help='Path to the file that will be analysed' )
parser.add_argument('--list', type=argparse.FileType('r'), required=True, help='File with a list of chros to be analysed (comma-separated)' )

args = parser.parse_args()

if __name__ == '__main__':
    
    #list of chros to analyse
    list_chros=[]
    for line in args.list:
        line=line.strip("\n")
        list_chros=line.split(',')
    
    bam = BamQC(args.filename,args.samtools)
    
    stats=bam.get_simple_stats()
    for k,v in stats.items(): print(k+":"+str(v))
    cov_list=bam.run_samtools_depth(chros=list_chros)
    aggr_cov=bam.aggregate_stats(cov_list)
    print("contig: %s" % aggr_cov.contig)
    print("breadth: %f" % aggr_cov.breadth)
    print("depth: %f" % aggr_cov.depth)
    print("max: %f" % aggr_cov.max)
    pass
