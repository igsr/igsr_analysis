import eHive
import os
import pdb
from VCF.VcfQC import VcfQC
from ReseqTrackDB import File
from ReseqTrackDB import ReseqTrackDB
from ReseqTrackDB import Attribute

class CollectVariantCallingMetrics(eHive.BaseRunnable):
    """run CollectVariantCallingMetrics on a VCF file"""
    
    def param_defaults(self):
        return {
        }

    def fetch_input(self):
        filepath = self.param_required('filepath')
        
        hostname=self.param_required('hostname')
        username=self.param_required('username')
        db=self.param_required('db')
        port=self.param_required('port')
        pwd=self.param_required('pwd')

        reseqdb = ReseqTrackDB(host=hostname,user=username,port=port,pwd=pwd,db=db)

        self.param('reseqdb', reseqdb)
         
    def run(self):
        
        reseqdb=self.param('reseqdb')
        
        self.warning('Analysing file: %s'% self.param_required('filepath'))

        vcf = VcfQC(vcf=self.param_required('filepath'),picard_folder=self.param_required('picard_folder'))

        cvcmetrics=""
        if self.param_is_defined('intervals'):
            cvcmetrics=vcf.run_CollectVariantCallingMetrics(outprefix=self.param_required('filepath'),truth_vcf=self.param_required('truth_vcf'),intervals=self.param('intervals'))
        else:
            cvcmetrics=vcf.run_CollectVariantCallingMetrics(outprefix=self.param_required('filepath'),truth_vcf=self.param_required('truth_vcf'))

        #store attributes
        if self.param_required('store_attributes')=='True':
            for attr,value in cvcmetrics.vc_detail_metrics.items():
                Attribute(table_name="file",other_id=fileO.dbID,name="CVCM_detail_"+attr,value=value).store(reseqdb)
        
            for attr,value in cvcmetrics.vc_summary_metrics.items():
                Attribute(table_name="file",other_id=fileO.dbID,name="CVCM_summary_"+attr,value=value).store(reseqdb)

        self.param('detail_metrics', cvcmetrics.vc_detail_metrics_file)
        self.param('summary_metrics', cvcmetrics.vc_summary_metrics_file)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'detail_metrics' : self.param('detail_metrics') }, 1)
        self.dataflow( { 'summary_metrics' : self.param('summary_metrics') }, 1)

