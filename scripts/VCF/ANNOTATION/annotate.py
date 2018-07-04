import argparse
import subprocess
import re
import os
import logging
import pdb
import re

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

def mod_AFs(AFfile):
    '''
    Function to process and select certain columns from the Allele Frequencies file calculated by function 'get_allele_frequencies'

    Args
    ----
    AFfile: string
            Path to file with allele frequencies calculated by function 'get_allele_frequencies'

    Returns
    -------
    New allele frequency file some modifications
    '''

    p = re.compile('0\.00+')
    AF_ixs=[7,10,13,16,19,22]

    pdb.set_trace()
    region=args.region.replace(':','.')
    pops=args.pops.replace(',','_')
    outfile="{0}/AFs.{1}.{2}.sub1.tsv".format(args.outdir,region,pops)
    wf=open(outfile,'w')
    with open(AFfile) as f:
        for line in f:
            line=line.rstrip("\n")
            b=line.split("\t")
            for i in AF_ixs:
                if p.match(b[i]):
                    b[i]='0'
            str="{0}-{1}\t{0}\t{1}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}".format(b[0],b[1],b[3],b[4],b[20],b[21],b[22],b[7],b[10],b[13],b[16],b[19])
            wf.write(str+"\n")
            
    wf.close()
    
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

    cmd="bcftools query -f '%CHROM-%POS\\t%INFO/{0}\\n' -r {1} {2} -o {3}".format(ann,args.region,args.ann_vcf,outfile)
    
    try:
        stdout = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        raise
        
    return outfile

def concat_tables(depth_f,ann_tab):
    '''
    Function to join 2 files by the common column (chr-pos in this case)

    Args
    ----
    depth_f: string
             Path to file with depths obtained after running 'get_annotation'
    ann_tab: string
             Path to file with annotations obtained by running the function 'mod_AFs'

    '''

    region=args.region.replace(':','.')
    tmpfile="{0}/concat_f.{1}.tmp.txt".format(args.outdir,region)
    cmd="join <(sort {0}) <(sort {1}) > {2}".format(depth_f,ann_tab,tmpfile)
    
    try:
        stdout = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        raise

    outfile="{0}/concat_f.{1}.txt".format(args.outdir,region)    
    wf=open(outfile,'w')

    wf.write("#CHR\tFROM\tTO\tREF\tALT\tDP\tAN\tAC\tAF\tEAS_AF\tEUR_AF\tAFR_AF\tAMR_AF\tSAS_AF\n")
    with open(tmpfile) as f:
        for line in f:
            line=line.rstrip("\n")
            b=line.split(" ")
            str="{0}\t{1}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}".format(b[2],b[3],b[5],b[6],b[1],b[7],b[8],b[9],b[10],b[11],b[12],b[13],b[14])
            wf.write(str+"\n")

    wf.close()

    os.remove(tmpfile)

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
logging.info("Running mod_AFs")
#modAFs=mod_AFs(AFs)
logging.info("Done!")
logging.info("Running get_annotation")
#depth_f=get_annotation(args.ann_vcf,'DP')
logging.info("Done!")
logging.info("Running concat_tables")
#concat_tab=concat_tables(depth_f,modAFs)
logging.info("Done!")
logging.info("Running get_overlapping_variants")
#annot_tab=get_overlapping_variants(concat_tab,args.exome,'EX_TARGET')
logging.info("Done!")
logging.info("Running add_annotation")
#annot_tab1=add_annotation(annot_tab,'VT','SNP')
logging.info("Done!")
annot_tab2=add_number_of_samples(annot_tab1,args.ann_vcf)

#delete old files
os.remove(AFs)
os.remove(modAFs)
os.remove(depth_f)
os.remove(concat_tab)
os.remove(annot_tab)
os.remove(annot_tab1)
