# Utility to run VerifyBAMID on a BAMs that is tracked in a ReseqTrackDB
# This script assumes that the filename contains the population ID for this particular sample BAM
# For example: NA20757.alt_bwamem_GRCh38DH_tags_stripped.20150718.TSI.low_coverage.chr20.bam contains
# the population id (TSI)

from BamQC import BamQC

from pprint import pprint
from ReseqTrackDB import *
import argparse
import os
import glob
import fnmatch
import re
import warnings
import logging

#get command line arguments

#RESEQTRACK DB conn params
parser = argparse.ArgumentParser(description='Script to calculate different QC metrics on BAMs')
parser.add_argument('--hostname', type=str, required=True, help='Hostname for ReseqTrack DB' )
parser.add_argument('--username', type=str, required=True, help='User for ReseqTrack DB' )
parser.add_argument('--port', type=int, required=True, help='Port number in the ReseqTrack DB' )
parser.add_argument('--pwd', type=str, help='PWD for the ReseqTrack DB' )
parser.add_argument('--db', type=str, required=True, help='DB name in the ReseqTrack DB' )

parser.add_argument('--exe', type=str, help='Folder containing the VerifyBamID binary' )
parser.add_argument('--filename', type=str, required=True, help='Base name of file that will be analysed' )
parser.add_argument('--outdir', type=str, required=True, help='output dir where results file will be written' )
parser.add_argument('--genotypes', type=str, required=True, help='folder with genotype files split by population' )
parser.add_argument('--prefix', type=str, required=True, help='prefix for attribute_name in the ReseqtrackDB attribute table' )

args = parser.parse_args()

if __name__ == '__main__':
    log_filename="%s/verifybam_id_%s.log" % (args.outdir,args.filename)

    logger = logging.getLogger("verifybamid")
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

    fileO=reseqdb.fetch_file_by_filename(args.filename)

    #basename for output
    basename=os.path.splitext(args.filename)[0]

    population_id=basename.split(".")[3]
    
    if not population_id:
        raise Exception("Population could not be extracted from %s:" % basename)

    genotype_f=glob.glob(args.genotypes+"/"+population_id+"*")

    if not len(genotype_f)==1:
        logger.warning("No population genotype file for %s" % basename)
        logger.info("Done!")
        exit(0)

    bam = BamQC(bam=fileO.path,verifybamid_folder=args.exe)

    outfiles=bam.run_verifybamid(genotype_f[0],basename,outdir=args.outdir)

    for file in outfiles:
        index=[]
        if fnmatch.fnmatch(file, '*.selfSM') or fnmatch.fnmatch(file, '*.selfRG'):
            with open(file) as f:
                for line in f:
                     line=line.rstrip('\n')
                     if re.match(r"^#", line):
                         bits=line.split('\t')
                         continue
                     else:
                         columns=line.split('\t')
                         sample=columns[bits.index('#SEQ_ID')]
                         rg=columns[bits.index('RG')]
                         freemix=columns[bits.index('FREEMIX')]
                         chipmix=columns[bits.index('CHIPMIX')]
                         name1="%s_%s_%s:freemix" %(args.prefix,sample,rg)
                         name2="%s_%s_%s:chipmix" %(args.prefix,sample,rg)
                         Attribute(table_name="file",other_id=fileO.dbID,name=name1,value=freemix).store(reseqdb)
                         Attribute(table_name="file",other_id=fileO.dbID,name=name2,value=chipmix).store(reseqdb)

    logger.info("Done!")
