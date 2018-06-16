'''
Utility to run all BAM QC tests as performed for the IGSR project over a BAM file/s stored in
a ReseqTrackDB

Created on 08 Feb 2017

@author: ernesto lowy: ernesto@ebi.ac.uk
'''

from BamQC import BamQC

from pprint import pprint
from ReseqTrackDB import *
import argparse
import os

#get command line arguments

#RESEQTRACK DB conn params
parser = argparse.ArgumentParser(description='Script to calculate different QC metrics on BAMs')
parser.add_argument('--hostname', type=str, required=True, help='Hostname for ReseqTrack DB' )
parser.add_argument('--username', type=str, required=True, help='User for ReseqTrack DB' )
parser.add_argument('--port', type=int, required=True, help='Port number in the ReseqTrack DB' )
parser.add_argument('--pwd', type=str, help='PWD for the ReseqTrack DB' )
parser.add_argument('--db', type=str, required=True, help='DB name in the ReseqTrack DB' )

parser.add_argument('--chkindel_exe', type=str, help='Folder containing the chk_indel_rg binary' )
parser.add_argument('--filename', type=str, required=True, help='Base name of file that will be analysed' )
parser.add_argument('--outdir', type=str, required=True, help='output dir where results file will be written' )

args = parser.parse_args()

def run_chkindel_rg (file):
    '''
    Function to run chkindel_rg over a BAM file.
    It will stored the results of the program in the Attribute table
    of the ReseqTrackDB.

    All attributes stored in the Attribute table will start with 'CHK_INDEL_' as the attribute_name

    Arguments
    ---------
    file: ReseqTrackDB file object

    Returns
    ------
    None

    '''

    bam = BamQC(bam=file.path,chk_indel_folder=args.chkindel_exe)

    #basename for output
    basename=os.path.splitext(args.filename)[0]

    listO=bam.run_chk_indel_rg(outfile=args.outdir+"/"+basename+".indel_rg.txt")

    for obj in listO:
        inser=obj.ins_in_short_homopolymer
        delet=obj.del_in_short

        #increment inser and delet counts by 1 to avoid divisions by 0
        if delet==0:
            delet+=1
            inser+=1
        #calculating ratio: ins-in-short-homopolymer/del-in-short and checking if it is greater than 5
        outcome="FAILED" if (float(inser)/float(delet))>5 else "PASS"
#        Attribute(table_name="file",other_id=file.dbID,name="CHK_INDEL_"+obj.RG+":ratio",value=outcome).store(reseqdb)
#        for attr, value in obj.__dict__.iteritems():
#            Attribute(table_name="file",other_id=file.dbID,name="CHK_INDEL_"+obj.RG+":"+attr,value=value).store(reseqdb)

if __name__ == '__main__':
    hostname=args.hostname
    username=args.username
    db=args.db
    port=args.port
    pwd=args.pwd

    reseqdb = ReseqTrackDB(host=hostname,user=username,port=port,pwd=pwd,db=db)

    file=reseqdb.fetch_file_by_filename(args.filename)

#    run_chkindel_rg(file)

