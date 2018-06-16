from ReseqTrackDB import *
import re
import pandas as pd
import numpy as np
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description='This script will analyze the different attributes in the Attribute table that have been calculated using the PyHive::PipeConfig::QC::RunBcfToolsStats pipeline over a certain VCF file. \
This pipeline needs be run with the "filter_str" parameter of the PyHive.VcfQC.BcftoolsStats runnable')

#RESEQTRACK DB conn params
parser.add_argument('--hostname', type=str, required=True, help='Hostname for ReseqTrack DB' )
parser.add_argument('--username', type=str, required=True, help='User for ReseqTrack DB' )
parser.add_argument('--port', type=int, required=True, help='Port number in the ReseqTrack DB' )
parser.add_argument('--pwd', type=str, help='PWD for the ReseqTrack DB' )
parser.add_argument('--db', type=str, required=True, help='DB name in the ReseqTrack DB' )

#files that will be analyzed
parser.add_argument('-vcfl','--vcflist', nargs='+', help='<Required> List of VCFs to analyze', required=True)
parser.add_argument('-cnames','--columnames', nargs='+', help='<Required> List of columns used in the .csv and .png files', required=True)
parser.add_argument('--chros_f', type=str, required=True, help='File with chros' )
parser.add_argument('--outprefix', type=str, required=True, help='Prefix for output files' )

args = parser.parse_args()

reseqtrackdb_object=ReseqTrackDB(host=args.hostname,user=args.username,port=args.port,pwd=args.pwd,db=args.db)

def get_annotation(attr_list,dict_str,pattern):

    patt = re.compile( pattern )
    patt1 = re.compile('^STATS_.*(chr.+)_filt' )
    for a in attr_list:
        match = patt.match( a.name )
        match1 = patt1.match( a.name )
        if match:
            if match1:
                dict_str[match1.group(1)].append(a.value)
            else:
                raise Exception("Error matching the attribute name")

    return dict_str

def get_chr_list(chros_f):

    chr_list=[]
    with open(chros_f) as f:
        for line in f:
            line=line.rstrip('\n')
            chr_list.append(line.split('\t')[0])
    return chr_list

def get_DF(chr_list,outprefix,thedict):
    data=[]
    for chrom in chr_list:
        if chrom in thedict:
            values=thedict[chrom]
            values=[chrom]+values
            data.append(values)
            
    DF=pd.DataFrame(data,columns=args.columnames).set_index('chrom')
    DF=DF.fillna(value=0)
    DF=DF.astype(float)
    ax=DF.plot(title=outprefix,rot=90,grid=True,figsize=(15,10))
    ax.set_xticks(range(len(DF.index)))
    ax.set_xticks(range(len(DF.index)))
    ax.set_xticklabels(DF.index)
    fig = ax.get_figure()
    fig.savefig("{0}.png".format(outprefix))
    DF.to_csv(outprefix+".csv",index=True,header=True)

chr_list=get_chr_list(args.chros_f)

dictA = defaultdict(list)
dictB = defaultdict(list)
dictC = defaultdict(list)
dictD = defaultdict(list)
dictE = defaultdict(list)

for f in args.vcflist:
    fileO=reseqtrackdb_object.fetch_file_by_url(url=f)
    attr_list=reseqtrackdb_object.fetch_attributes_by_other_id(fileO.dbID)
    dictA=get_annotation(attr_list=attr_list,dict_str=dictA,pattern='^STATS_chr.+_filt_number of records:')
    dictB=get_annotation(attr_list=attr_list,dict_str=dictB,pattern='^STATS_chr.+_filt_number of SNPs:')
    dictC=get_annotation(attr_list=attr_list,dict_str=dictC,pattern='^STATS_chr.+_filt_number of indels:')
    dictD=get_annotation(attr_list=attr_list,dict_str=dictD,pattern='^STATS_chr.+_filt_number of multiallelic sites:')
    dictE=get_annotation(attr_list=attr_list,dict_str=dictE,pattern='^STATS_ts_tv_chr.+_filt')

get_DF(chr_list=chr_list,outprefix='{0}_no_of_records'.format(args.outprefix),thedict=dictA)
get_DF(chr_list=chr_list,outprefix='{0}_no_of_SNPs'.format(args.outprefix),thedict=dictB)
get_DF(chr_list=chr_list,outprefix='{0}_no_of_INDELs'.format(args.outprefix),thedict=dictC)
get_DF(chr_list=chr_list,outprefix='{0}_no_of_multiallelic_sites'.format(args.outprefix),thedict=dictD)
get_DF(chr_list=chr_list,outprefix='{0}_ts_tv'.format(args.outprefix),thedict=dictE)
