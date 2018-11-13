'''
Created on 25 Oct 2016

@author: ernesto
'''

from ReseqTrackDB import ReseqTrackDB
from ReseqTrackDB import Attribute
from BamQC import BamQC

if __name__ == '__main__':
    hostname = None
    username = None
    db = None
    port = None
    pwd = None

    reseqdb = ReseqTrackDB(host=hostname, user=username, port=port, pwd=pwd, db=db)

    file = reseqdb.fetch_file_by_filename("HG03078.alt_bwamem_GRCh38DH_"
                                          "tags_stripped.20150826.MSL.exome.chr20.bam")

    #check if get_simple_stats has been already run on this file
    attrb = reseqdb.fetch_attribute_by_params('file', file.dbID, 'no_properly_paired')

    if attrb:
        print("Exists!")

    bam = BamQC(bam=file.path, samtools_folder="/Users/ernesto/lib/samtools-1.3.1/")

    stats = bam.get_simple_stats()
    for k, v in stats.items():
        Attribute(table_name="file", other_id=file.dbID, name=k, value=v).store(reseqdb)
