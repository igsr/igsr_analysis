import eHive
import os
from VcfFilter import VcfFilter

class SplitVariants(eHive.BaseRunnable):
    """Split VCF into SNPs and INDELs"""

    def run(self):

        self.warning('Analysing file: %s'% self.param_required('filepath'))
        
        filepath=self.param_required('filepath')

        vcf = VcfFilter(vcf=filepath, bcftools_folder=self.param_required('bcftools_folder'))

        outprefix=""

        if self.param_is_defined('work_dir'):
            file=os.path.split(filepath)[1]
            outprefix=self.param('work_dir')+"/"+file
        else:
            outprefix=filepath

        outfile=vcf.filter_by_variant_type(type=self.param_required('type'),outprefix=outprefix)
        
        self.param('out_vcf', outfile)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'out_vcf' : self.param('out_vcf') }, 1)

