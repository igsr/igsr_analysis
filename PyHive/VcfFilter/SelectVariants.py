import eHive
import os
import pdb
from VCF.VCFfilter.BCFTools import BCFTools

class SelectVariants(eHive.BaseRunnable):
    """Select the variants from a VCF file"""

    def run(self):

        self.warning('Analysing file: %s'% self.param_required('filepath'))

        filepath = self.param_required('filepath')

        vcf = BCFTools(vcf=filepath, bcftools_folder=self.param_required('bcftools_folder'))

        outprefix = ""
        if self.param_is_defined('work_dir'):
            file = os.path.split(filepath)[1]
            outprefix = self.param('work_dir')+"/"+file
        else:
            outprefix = filepath

        uncalled = None
        if self.param_is_defined('uncalled'):
            uncalled = self.param('uncalled')
            if uncalled != 'include' and uncalled != 'exclude':
                raise Exception("Valid 'uncalled' values are 'include'/'exclude'")

        threads = 1
        if self.param_is_defined('threads'):
            threads = self.param('threads')

        outfile = vcf.select_variants(outprefix=outprefix,
                                      threads=threads,
                                      uncalled=uncalled)
        
        self.param('out_vcf', outfile)

    def write_output(self):
        self.warning('Work is done!')

        self.dataflow({'out_vcf': self.param('out_vcf')}, 1)
