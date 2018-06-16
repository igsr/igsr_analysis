import eHive
import os
import pdb
from VCF.VCFfilter.BCFTools import BCFTools

class SplitVariants(eHive.BaseRunnable):
    """Split VCF into SNPs and INDELs"""

    def run(self):

        self.warning('Analysing file: %s'% self.param_required('filepath'))
        
        filepath=self.param_required('filepath')

        vcf = BCFTools(vcf=filepath, bcftools_folder=self.param_required('bcftools_folder'))

        outprefix=""

        if self.param_is_defined('work_dir'):
            file=os.path.split(filepath)[1]
            outprefix=self.param('work_dir')+"/"+file
        else:
            outprefix=filepath

        biallelic=None
        if self.param_is_defined('biallelic'):
            if self.param('biallelic')=="True":
                biallelic=True
            elif self.param('biallelic')=="False":
                biallelic=False
            else:
                raise Exception("Error. biallelic option should be True or False")

        compress=None
        if self.param_is_defined('compress'):
            if self.param('compress')=="True":
                compress=True
            elif self.param('compress')=="False":
                compress=False
            else:
                raise Exception("Error. compress option should be True or False")

        outfile=vcf.filter_by_variant_type(v_type=self.param_required('type'),outprefix=outprefix, biallelic=biallelic, compress=compress)
        
        self.param('out_vcf', outfile)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'out_vcf' : self.param('out_vcf') }, 1)

