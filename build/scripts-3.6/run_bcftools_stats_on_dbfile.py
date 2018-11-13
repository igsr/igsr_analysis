from VcfQC import VcfQC
from ReseqTrackDB import File
from ReseqTrackDB import ReseqTrackDB
from ReseqTrackDB import Attribute

import argparse
import os
import logging

#get command line arguments

parser = argparse.ArgumentParser(description='Script to run bcftools stats on a certain VCF file')

'''
Reseqtrack DB connection parameters
'''
parser.add_argument('--hostname', type=str, required=True, help='Hostname for ReseqTrack DB' )
parser.add_argument('--username', type=str, required=True, help='User for ReseqTrack DB' )
parser.add_argument('--port', type=int, required=True, help='Port number in the ReseqTrack DB' )
parser.add_argument('--pwd', type=str, help='PWD for the ReseqTrack DB' )
parser.add_argument('--db', type=str, required=True, help='DB name in the ReseqTrack DB')

parser.add_argument('--bcftools_folder', type=str, required=True, help='Folder containing the VCFtools binary' )
parser.add_argument('--filename', type=str, required=True, help='Base name of the VCF file that will be analysed' )
parser.add_argument('--outdir', type=str, required=True, help='Directory used to put the output files.')

args = parser.parse_args()

if __name__ == '__main__':
    if os.path.isdir(args.outdir) == False:
        raise Exception("Output dir does not exist: %s"%args.outdir)
    log_filename="bcftools_stats_%s.log"% args.filename
    
    logger = logging.getLogger("bcftools_stats")
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
    l=os.path.basename(file.path).split('.')
    l.remove('vcf')
    l.remove('gz')
    outprefix='.'.join(l)
    
    vcfQC = VcfQC(vcf=file.path,bcftools_folder=args.bcftools_folder)
    stats=vcfQC.stats(outprefix=outprefix,outdir=args.outdir)

    #store attributes
    for attr,value in stats.summary_numbers.items():
        Attribute(table_name="file",other_id=file.dbID,name="STATS_"+attr,value=value).store(reseqdb)
    
    Attribute(table_name="file",other_id=file.dbID,name="STATS_ts_tv",value=stats.ts_tv).store(reseqdb)
    Attribute(table_name="file",other_id=file.dbID,name="STATS_ts_tv_1stalt",value=stats.ts_tv_1stalt).store(reseqdb)
    Attribute(table_name="file",other_id=file.dbID,name="STATS_no_singleton_snps",value=stats.no_singleton_snps).store(reseqdb)

    #store file
    stats_f=File(path=stats.filename,type=file.type+"_STATS",host_id=1,withdrawn=0)
    stats_f.store(reseqdb,do_md5=True)
            
    logger.info("Done!.")
