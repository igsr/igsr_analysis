import eHive
import os

from VCFIntegration.SNPTools import SNPTools

class SNPTools_bamodel(eHive.BaseRunnable):
    """Run SNPTools bamodel on a VCF containing biallelic SNPs"""

    def run(self):
        vcf_g = SNPTools(vcf=self.param_required('vcf_file'),
                         snptools_folder=self.param_required('snptools_folder'))

        verbose = None
        if self.param_is_defined('verbose'):
            verbose = True
        else:
            verbose = False

        if self.param_is_defined('work_dir'):
            if not os.path.isdir(self.param('work_dir')):
                os.makedirs(self.param('work_dir'))

        raw_f = vcf_g.run_bamodel(sample=self.param_required('sample'),
                                  bamfiles=self.param_required('bamlist'),
                                  outdir=self.param_required('work_dir'),
                                  verbose=verbose)

        self.param('raw_f', raw_f)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow({'raw_f': self.param('raw_f')}, 1)
