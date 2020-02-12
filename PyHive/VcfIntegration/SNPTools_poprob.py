import eHive
import os
from VCFIntegration.SNPTools import SNPTools

class SNPTools_poprob(eHive.BaseRunnable):
    """Run SNPTools poprob on a VCF containing biallelic SNPs"""

    def run(self):
        vcf_g = SNPTools(vcf=self.param_required('vcf_file'),
                         snptools_folder=self.param_required('snptools_folder'))

        outprefix = os.path.split(self.param_required('outprefix'))[1]

        if self.param_is_defined('work_dir'):
            if not os.path.isdir(self.param('work_dir')):
                os.makedirs(self.param('work_dir'))

        prob_f = ""
        if self.param_is_defined('verbose'):
            prob_f = vcf_g.run_poprob(outprefix=outprefix,
                                      rawlist=self.param_required('rawlist'),
                                      outdir=self.param_required('work_dir'),
                                      verbose=True)
        else:
            prob_f = vcf_g.run_poprob(outprefix=outprefix,
                                      rawlist=self.param_required('rawlist'),
                                      outdir=self.param_required('work_dir'),
                                      verbose=False)

        self.param('prob_f', prob_f)

    def write_output(self):
        self.warning('Work is done!')

        self.dataflow({'prob_f': self.param('prob_f')}, 1)
