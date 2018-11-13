from BamQC import BamQC,CMetrics

from pprint import pprint
from ReseqTrackDB import *
import argparse
import os

#get command line arguments

#RESEQTRACK DB conn params
parser = argparse.ArgumentParser(description='Script to calculate different QC metrics on the WGS BAMs')
parser.add_argument('--hostname', type=str, required=True, help='Hostname for ReseqTrack DB' )
parser.add_argument('--username', type=str, required=True, help='User for ReseqTrack DB' )
parser.add_argument('--port', type=int, required=True, help='Port number in the ReseqTrack DB' )
parser.add_argument('--pwd', type=str, help='PWD for the ReseqTrack DB' )
parser.add_argument('--db', type=str, required=True, help='DB name in the ReseqTrack DB' )

parser.add_argument('--picard', type=str, required=False, help='Folder containing the Picard jar file' )
parser.add_argument('--java', type=str, required=False, help='Folder containing the Java binary' )
parser.add_argument('--samtools', type=str, required=True, help='Folder containing the SAMtools binary' )
parser.add_argument('--filename', type=str, required=True, help='Base name of file that will be analysed' )
parser.add_argument('--outdir', type=str, required=True, help='output dir where metrics and coverage barplot files will be created' )
parser.add_argument('--targetfile', type=str, required=True, help='File containing exome coordinates' )

args = parser.parse_args()

if __name__ == '__main__':
    hostname=args.hostname
    username=args.username
    db=args.db
    port=args.port
    pwd=args.pwd
    
    reseqdb = ReseqTrackDB(host=hostname,user=username,port=port,pwd=pwd,db=db)
    
    file=reseqdb.fetch_file_by_filename(args.filename)
    
    bam = BamQC(bam=file.path,samtools_folder=args.samtools,java_folder=args.java,picard_folder=args.picard)
    
    list_of_samples=bam.list_of_samples()
    list_of_readgroups=bam.list_of_readgroups()
    stats=bam.get_simple_stats()
    for k,v in stats.items(): Attribute(table_name="file",other_id=file.dbID,name=k,value=v).store(reseqdb)
    for i in range(len(list_of_readgroups)): Attribute(table_name="file",other_id=file.dbID,name="readgroup%d" % (i+1),value=list_of_readgroups[i]).store(reseqdb)
    for i in range(len(list_of_samples)): Attribute(table_name="file",other_id=file.dbID,name="sample%d" % (i+1),value=list_of_samples[i]).store(reseqdb)
    
    #basename for output
    basename=os.path.splitext(args.filename)[0]
    CHsMetricsO=bam.run_CollectHsMetrics(args.targetfile,outfile=args.outdir+"/"+basename+".CollectHsMetrics.metrics.txt")
    for k,v in CHsMetricsO.metrics.items(): Attribute(table_name="file",other_id=file.dbID,name="PICARD:"+k,value=v).store(reseqdb)
    #edit xlim and ylim depending on coverage levels
    CHsMetricsO.create_cov_barplot(args.outdir+"/"+basename+".CollectHsMetrics.barplot.pdf")
    
