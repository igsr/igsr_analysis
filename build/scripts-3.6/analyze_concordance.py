from VcfQC import VcfQC
from VcfQC import GTPconcordance
from ReseqTrackDB import ReseqTrackDB
from ReseqTrackDB import Attribute

import argparse
import os
import logging

#get command line arguments

parser = argparse.ArgumentParser(description='Script to run assess the genotype concordance between 2 VCFs, one of them considered as the truth set')

'''
Reseqtrack DB connection parameters
'''
parser.add_argument('--hostname', type=str, required=True, help='Hostname for ReseqTrack DB' )
parser.add_argument('--username', type=str, required=True, help='User for ReseqTrack DB' )
parser.add_argument('--port', type=int, required=True, help='Port number in the ReseqTrack DB' )
parser.add_argument('--pwd', type=str, help='PWD for the ReseqTrack DB' )
parser.add_argument('--db', type=str, required=True, help='DB name in the ReseqTrack DB')

parser.add_argument('--working_dir', type=str, required=True, help='Working dir' )
parser.add_argument('--picard_folder', type=str, required=True, help='Picard jar file' )
parser.add_argument('--filename', type=str, required=True, help='VCF file in the DB with call set' )
parser.add_argument('--truth_vcf', type=str, required=True, help='VCF file with truth set' )
parser.add_argument('--sample_f', type=str, required=True, help='File containing the samples to be analysed' )
parser.add_argument('--prefix', type=str, required=True, help='Prefix for output files' )
parser.add_argument('--intervals', type=str, required=False, help='Interval file' )


args = parser.parse_args()

if __name__ == '__main__':
    if os.path.isdir(args.working_dir) == False:
        raise Exception("Dir does not exist: %s"%args.working_dir)

    #getting basename for --sample_f
    basename=os.path.basename(args.sample_f).split('.')[0]
    log_filename=args.working_dir+"/genotype_concordance_%s.log"% basename
    
    logger = logging.getLogger("gt_concordance")
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

    file=reseqdb.fetch_file_by_filename(args.filename)

    vcfQC = VcfQC(vcf=file.path,picard_folder=args.picard_folder)

    if os.path.isfile(args.sample_f) == False:
        raise Exception("File does not exist")

    with open(args.sample_f) as f:
        for s in f:
            s=s.rstrip('\n')
            gtp_con=vcfQC.calc_concordance(truth_vcf=args.truth_vcf,truth_sample=s,call_sample=s,outprefix=s+"_"+args.prefix, intervals=args.intervals)

            #store attributes
            for attr,value in gtp_con.summary_metrics_snps.items():
                Attribute(table_name="file",other_id=file.dbID,name="GT_CONC_"+s+"_"+attr,value=value).store(reseqdb)
                
    logger.info("Done!.")
