import argparse
import subprocess
import re
import os
import logging
import pdb
import re
import pandas as pd


parser = argparse.ArgumentParser(description='Script to annotate a VCF')

parser.add_argument('--AFcalc',help='Folder containing the Perl script (calculate_allele_frq_from_vcf.pl) used to calculate the allele frequency',required=True)
parser.add_argument('--phased_vcf',help='Path to phased VCF that will be annotated',required=True)
parser.add_argument('--sample_panel',help='Panel file with population/superpopulation for each sample',required=True)
parser.add_argument('--tabix',help='Folter containing the tabix binary',required=True)
parser.add_argument('--region',help='Region to analyse. i.e. chr20:1-1000',required=True)
parser.add_argument('--pops',help='Populations or superpopulations to analyse. i.e. EAS,EUR,AFR,AMR,SAS',required=True)
parser.add_argument('--exome',help='BED file with the exome pull down targets',required=True)
parser.add_argument('--outdir',help='Folder used to put the output files',required=True)
parser.add_argument('--ann_vcf',help='Path to VCF with annotations (i.e. DEPTH) that will be used to create the final annotation table',required=True)

args = parser.parse_args()

logging.basicConfig(filename="sample.log", filemode="w", level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')

def get_allele_frequencies():
    '''
    Function to calculate AFs on a VCF file

    Returns
    -------
    Path to a tsv file containing the allele frequencies

    '''

    region=args.region.replace(':','.')
    pops=args.pops.replace(',','_')
    outfile="{0}/AFs.{1}.{2}.tsv".format(args.outdir,region,pops)

    cmd="perl {0}/calculate_allele_frq_from_vcf.pl -vcf {1} -sample_panel {2} -out_file {3} -region {4} -tabix {5} -pop {6}".format(args.AFcalc,
                                                                                                                                    args.phased_vcf,
                                                                                                                                    args.sample_panel,
                                                                                                                                    outfile,
                                                                                                                                    args.region,
                                                                                                                                    args.tabix,
                                                                                                                                    args.pops)
    try:
        stdout = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        raise
    
    return outfile

def get_annotation(vcf, ann):
    '''
    Function to run bcftools query on a VCF file to retrieve a particular annotation

    Args
    ----
    vcf: string
         Path to VCF file containing the annotations (i.e. DP) that will included in the annotations table
    ann: string
         Annotation that will be extracted from vcf (i.e. DP)

    Returns
    -------
    A file containing the following 2 columns:
    <chr-pos>\t<DP>

    '''

    region=args.region.replace(':','.')
    outfile="{0}/depths.{1}.txt".format(args.outdir,region)

    cmd="bcftools query -f '%CHROM\\t%POS\\t%INFO/{0}\\n' -r {1} {2} -o {3}".format(ann,args.region,args.ann_vcf,outfile)

    
    try:
        stdout = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        raise
        
    return outfile

def process_AFs(AF_f,depth_f):
    '''
    Function to create a table containing the AFs and the depths for each position

    Args
    ----
    AF_f: string
          Path to file containing the allele frequencies created by 'get_allele_frequencies'
    depth_f: string
             Path to file containing the depths per position created by function 'get_annotation'

    Returns
    -------
    Path to table containing the AFs and the depths
    '''
 
    #Process the AFs
    cols=['#CHROM','POS','REF','ALT','ALL_TOTAL_CNT','ALL_ALT_CNT','ALL_FRQ','EAS_FRQ','EUR_FRQ','AFR_FRQ','AMR_FRQ','SAS_FRQ']
    names=['#CHR','POS', 'REF','ALT','AN','AC','AF','EAS_AF','EUR_AF','AFR_AF','AMR_AF','SAS_AF']
    final=['#CHR','FROM','TO','REF','ALT','DP','AN','AC','AF','EAS_AF','EUR_AF','AFR_AF','AMR_AF','SAS_AF']
    AFs_df=pd.read_csv(AF_f, sep="\t", usecols=cols, index_col=False)[cols]
    AFs_df.columns =names
    AFs_df=AFs_df.round(decimals=2).astype(object)

    #Proces the depths
    depths_DF=pd.read_csv(depth_f, sep="\t",  header=None)
    depths_DF.columns=["#CHR","POS","DP"]

    #Merge by chr-pos
    merged=pd.merge(AFs_df, depths_DF, on=['#CHR','POS'], how='inner')
    mergedA=merged.assign(TO = lambda x: x.POS)
    mergedA=mergedA.rename(columns={'POS': 'FROM'})
    mergedFinal=mergedA[final]

    #output filename
    region=args.region.replace(':','.')
    pops=args.pops.replace(',','_')

    outfile="{0}/modAF.{1}.{2}.txt".format(args.outdir,pops,region)
    mergedFinal.to_csv(outfile,sep='\t',index=False,float_format='%g')

    return outfile
    

def get_overlapping_variants(file1,file2,label):
    '''
    Function to identify the intersecting variants between file1 and file2. This function makes use of Bedtools intersect

    Args
    ----
    file1: string
           Path to 1st file containing the variants to be intersected. A new column will be appended to this file depending on the intersection 
           with the 2nd file
    file2: string
           Path to 2nd file containing the variants to be intersected
    label: string
           Label used for the column added to file1 when 2 variants intersect

    Returns
    -------
    File1 will be returned with a column added depending on the intersection with file1
    '''

    cmd="""awk '{ print $1\"\\t\"$2\"\\t\"$3 }' %s | bedtools intersect -a - -b %s """ % (file1,file2)

    intersect_list=[]

    try:
        out = subprocess.check_output(cmd, shell=True,executable='/bin/bash').decode("utf-8").split("\n")
        for l in out:
            if len(l.split('\t'))==3: intersect_list.append(l.split('\t')[1])
    except subprocess.CalledProcessError as e:
        raise

    region=args.region.replace(':','.')
    outfile="{0}/annot_tab.{1}.txt".format(args.outdir,region)
    wf=open(outfile,'w')
    with open(file1) as f:
        for line in f:
            line=line.rstrip("\n")
            if line.startswith("#"): 
                line=line+"\t{0}".format(label)
                wf.write(line+"\n")
                continue
            if line.split("\t")[1] in intersect_list:
                line= line+"\t1"
            else:
                line=line+"\t0"
            wf.write(line+"\n")
    wf.close()
    
    return outfile

def add_annotation(file1,header,label):
    '''
    Function to add a particular label to all rows in the table

    Args
    ----
    file1: string
           Path to table that will be used in order to add the annotation as the last column
    header: string
            Name used in the header for this annotation
    label: string
           Label that will be used and will be added as the last column in the table

    Returns
    -------
    Path to the file that will contain the additional column
    '''

    region=args.region.replace(':','.')
    outfile="{0}/annot_tab1.{1}.txt".format(args.outdir,region)

    wf=open(outfile,'w')
    with open(file1) as f:
        for line in f:
            line=line.rstrip("\n")
            if line.startswith('#'):
                line=line+"\t{0}".format(header)
                wf.write(line+"\n")
                continue
            line=line.rstrip("\n")
            line= line+"\t{0}".format(label)
            wf.write(line+"\n")
    wf.close()

    return outfile

def add_number_of_samples(file1,ann_vcf):
    '''
    Function to add the NS (number of samples) annotation to the annotation table

    Args
    ----
    file1: string
           Path to table that will be used in order to add the annotation as the last column
    ann_vcf: string
             Path to VCF with annotations that will be used to obtain the number of samples
   
    Returns
    -------
    Path to the file that will contain the additional column
             
    '''
    
    region=args.region.replace(':','.')
    outfile="{0}/annot_tab2.{1}.txt".format(args.outdir,region)


    cmd="bcftools query -l {0} |wc -l".format(ann_vcf)

    try:
        number_of_samples = subprocess.check_output(cmd, shell=True).decode("utf-8").rstrip('\n')
        wf=open(outfile,'w')
        with open(file1) as f:
            for line in f:
                line=line.rstrip("\n")
                if line.startswith('#'):
                    line=line+"\tNS"
                    wf.write(line+"\n")
                    continue
                else:
                    line= line+"\t{0}".format(number_of_samples)
                    wf.write(line+"\n")
    except subprocess.CalledProcessError as e:
        raise

logging.info("Starting the script")
logging.info("Running get_allele_frequencies")
AFs=get_allele_frequencies()
logging.info("Done!")
logging.info("Running get_annotation")
depth_f=get_annotation(args.ann_vcf,'DP')
logging.info("Done!")
logging.info("Running process_AFs")
afs_table=process_AFs(AFs,depth_f)
logging.info("Running get_overlapping_variants")
annot_tab=get_overlapping_variants(afs_table,args.exome,'EX_TARGET')
logging.info("Done!")
logging.info("Running add_annotation")
annot_tab1=add_annotation(annot_tab,'VT','SNP')
logging.info("Done!")
logging.info("Running add_number_of_samples")
annot_tab2=add_number_of_samples(annot_tab1,args.ann_vcf)
logging.info("Done!")

#delete old files
#os.remove(AFs)
#os.remove(afs_table)
#os.remove(depth_f)
#os.remove(annot_tab)
#os.remove(annot_tab1)
