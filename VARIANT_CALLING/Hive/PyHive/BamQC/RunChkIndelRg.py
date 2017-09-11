import eHive
import os
from BamQC import BamQC
from ReseqTrackDB import *

class RunChkIndelRg(eHive.BaseRunnable):
    """run chk_indel_rg on a BAM file"""
    
    def param_defaults(self):
        return {
        }

    def fetch_input(self):
        hostname=self.param('hostname')
        username=self.param('username')
        db=self.param('db')
        port=self.param('port')
        pwd=self.param('pwd')

        reseqdb = ReseqTrackDB(host=hostname,user=username,port=port,pwd=pwd,db=db)

        self.param('reseqdb', reseqdb)

    def run(self):
        filepath=self.param_required('filepath')
        reseqdb=self.param('reseqdb')

        file=reseqdb.fetch_file_by_url(filepath)

        self.warning('Analysing file: %s'% filepath)
        
        bam = BamQC(bam=filepath,chk_indel_folder=self.param_required('chk_indel_rg_folder'))

        outfile=None

        if self.param_is_defined('work_dir'):
            #basename for output
            basename=os.path.basename(filepath)
            outfile="{0}/{1}.indel_rg.txt".format(self.param('work_dir'), basename)
        else:
            (path,filename)=os.path.split(filepath)
            outfile="{0}/{1}.indel_rg.txt".format(path, filename)

        listO=bam.run_chk_indel_rg(outfile=outfile)

        chkindel_list=[]

        for obj in listO:
            outcome=obj.calc_ratio()
            for attr, value in obj.__dict__.items():
                attrb=Attribute(table_name="file",other_id=file.dbID,name="CHK_INDEL_"+obj.RG+":"+attr,value=value)
                chkindel_list.append(attrb.__dict__)
                            
        self.param('chkindel_list', chkindel_list)
        
    def write_output(self):
        self.warning('Work is done!')
        self.warning('{0} different Attributes were passed down'.format(len(self.param('chkindel_list'))))
        for attrb in self.param('chkindel_list'):
            self.dataflow( { 'attrb' : attrb }, 1)
