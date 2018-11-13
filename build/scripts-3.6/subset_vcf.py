
from VcfQC import VcfQC
from ReseqTrackDB import File
from ReseqTrackDB import ReseqTrackDB

import argparse
import os
import logging
import datetime

#get command line arguments

parser = argparse.ArgumentParser(description='Script to subset a VCF by excluding the variants within the regions defined by a BED file')

'''
Reseqtrack DB connection parameters
'''
parser.add_argument('--hostname', type=str, required=True, help='Hostname for ReseqTrack DB' )
parser.add_argument('--username', type=str, required=True, help='User for ReseqTrack DB' )
parser.add_argument('--port', type=int, required=True, help='Port number in the ReseqTrack DB' )
parser.add_argument('--pwd', type=str, help='PWD for the ReseqTrack DB' )
parser.add_argument('--db', type=str, required=True, help='DB name in the ReseqTrack DB' )
parser.add_argument('--type', type=str, required=True, help='Type of the new VCF file' )

parser.add_argument('--vcftools_folder', type=str, required=True, help='Folder containing the VCFtools binary' )
parser.add_argument('--bgzip_folder', type=str, required=True, help='Folder containing the bgzip binary')
parser.add_argument('--filename', type=str, required=True, help='Name (without the fullpath) of the VCF file that will be analysed. It assumes that the filename format is for example lc_bams.gatk.xxxx.vcf.gz, where lc_bams is the analysis group and gatk is the method used' )
parser.add_argument('--bed', type=str, required=True, help='BED file containing the coordinates to exclude' )
parser.add_argument('--outsuffix', type=str, required=True, help='Suffix for vcf output file. i.e. no_cms or no_offtarget' )
parser.add_argument('--outdir', type=str, required=True, help='Directory used to put the output files.' )

args = parser.parse_args()

if __name__ == '__main__':
    if os.path.isdir(args.outdir) == False:
        raise Exception("Output dir does not exist: %s"%args.outdir)
    
    hostname=args.hostname
    username=args.username
    db=args.db
    port=args.port
    pwd=args.pwd

    reseqdb = ReseqTrackDB(host=hostname,user=username,port=port,pwd=pwd,db=db)

    file=reseqdb.fetch_file_by_filename(args.filename)
    
    #constructing the out filename
    now = datetime.datetime.now().strftime('%Y%m%d')
    bits= os.path.basename(file.name).split('.')
    outprefix=bits[0]+"."+bits[1]+"."+args.outsuffix+"."+now

    log_filename="subset_vcf_%s.log"% outprefix

    logger = logging.getLogger("subset_vcf")
    logger.setLevel(logging.INFO)

    #  create the logging file handler
    fh = logging.FileHandler(log_filename)

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)

    # add handler to logger object
    logger.addHandler(fh)

    logger.info("Program started")

    vcfQC = VcfQC(vcf=file.path,bgzip_folder=args.bgzip_folder,vcftools_folder=args.vcftools_folder)
    vcffile=vcfQC.subset_vcf(bed=args.bed,outprefix=outprefix,outdir=args.outdir,create_index=True)

    f=File(path=vcffile,type=args.type,host_id=1,withdrawn=0)
    f.store(reseqdb,do_md5=True)
    
    logger.info("Done!.")
