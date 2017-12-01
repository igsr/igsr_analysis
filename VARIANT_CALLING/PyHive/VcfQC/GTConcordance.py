import eHive
import os
import datetime
from VcfQC import VcfQC
from ReseqTrackDB import File
from ReseqTrackDB import ReseqTrackDB
from ReseqTrackDB import Attribute

class GTConcordance(eHive.BaseRunnable):
    """Run Picard's GenotypeConcordance tool on a certain multi-sample VCF file"""

    def param_defaults(self):
        return {
        }

    def fetch_input(self):
        filename = self.param_required('filename')

        hostname=self.param_required('hostname')
        username=self.param_required('username')
        db=self.param_required('db')
        port=self.param_required('port')
        pwd=self.param_required('pwd')

        reseqdb = ReseqTrackDB(host=hostname,user=username,port=port,pwd=pwd,db=db)

        fileO=reseqdb.fetch_file_by_filename(filename)
        self.param('file_object', fileO)
        self.param('reseqdb', reseqdb)

    def run(self):
        fileO=self.param('file_object')
        reseqdb=self.param('reseqdb')
        sample=self.param('sample')

        self.warning('Analysing file: %s with sample: %s'% (fileO.name,sample))
        
        vcfQC = VcfQC(vcf=fileO.path,picard_folder=self.param_required('picard_folder'))
        gtp_con=vcfQC.calc_concordance(truth_vcf=self.param_required('truth_vcf'),truth_sample=sample,call_sample=sample,outprefix=sample+"_"+self.param_required('outprefix'),outdir=self.param_required('final_dir'),intervals=self.param('intervals'))
         
        if self.param_is_defined('store_attributes'):
            #store attributes
            for attr,value in gtp_con.summary_metrics_snps.items():
                Attribute(table_name="file",other_id=fileO.dbID,name="GT_CONC_"+sample+"_"+self.param_required('outprefix')+"_"+attr,value=value).store(reseqdb)

    def write_output(self):
        self.warning('Work is done!')

