import eHive
import os
from VCFfilter.BCFTools import BCFTools

class BcftoolsFilter(eHive.BaseRunnable):
    """run Bcftools filter on a VCF"""

    def run(self):

        self.warning('Analysing file: %s'% self.param_required('filepath'))
        
        vcf = BCFTools(vcf=self.param_required('filepath'),bcftools_folder=self.param_required('bcftools_folder'))

        outvcf=vcf.filter(name=self.param_required('filter_name'),expression=self.param_required('filter_expression'))
        
        self.param('out_vcf', outvcf)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'out_vcf' : self.param('out_vcf') }, 1)


