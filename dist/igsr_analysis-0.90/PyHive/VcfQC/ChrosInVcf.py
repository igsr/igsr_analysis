import eHive
import os
from VcfQC import VcfQC
from ReseqTrackDB import File
from ReseqTrackDB import ReseqTrackDB
from ReseqTrackDB import Attribute

class ChrosInVcf(eHive.BaseRunnable):
    """This runnable will compare the chros passed using the 'chr_file' option with the chros present in the VCF"""
    
    def param_defaults(self):
        return {
        }
    
    def run(self):

        filepath=self.param_required('filepath')

        self.warning('Analysing file: %s'% filepath)

        vcfQC = VcfQC(vcf=filepath,bcftools_folder=self.param_required('bcftools_folder'))
        chr_dict=vcfQC.get_chros(filter_str=self.param_required('filter_str'),chr_f=self.param_required('chr_file'))

        outcome="PASS"
        if len(chr_dict['in_A'])>0:
            outcome="FAIL"
            self.warning("TEST FAILED! There are chros present in {0} and not present in {1}".format(self.param_required('filepath'),self.param_required('chr_file')))
            self.warning("Chros are: {0}".format(','.join(chr_dict['in_A'])))

        if len(chr_dict['in_B'])>0:
            outcome="FAIL"
            self.warning("TEST FAILED! There are chros not present in {0} and present in {1}".format(self.param_required('filepath'),self.param_required('chr_file')))
            self.warning("Chros are: {0}".format(','.join(chr_dict['in_B'])))

        if self.param_required('store_attributes')=='True':
            hostname=self.param('hostname')
            username=self.param('username')
            db=self.param('db')
            port=self.param('port')
            pwd=self.param('pwd')

            reseqdb = ReseqTrackDB(host=hostname,user=username,port=port,pwd=pwd,db=db)

            fileO=reseqdb.fetch_file_by_url(filepath)
            
            #store attributes
            Attribute(table_name="file",other_id=fileO.dbID,name="QC_CHRCHECK",value=outcome).store(reseqdb)
            if len(chr_dict['in_A'])>0:
                Attribute(table_name="file",other_id=fileO.dbID,name="QC_CHRCHECK_IN_A",value=','.join(chr_dict['in_A'])).store(reseqdb)
            if len(chr_dict['in_B'])>0:
                Attribute(table_name="file",other_id=fileO.dbID,name="QC_CHRCHECK_IN_B",value=','.join(chr_dict['in_B'])).store(reseqdb)

        self.param('outcome', outcome)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'outcome' : self.param('outcome') }, 1)
