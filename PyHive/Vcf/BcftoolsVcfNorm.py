import eHive
import os
from VCF.VcfNormalize import VcfNormalize

class BcftoolsVcfNorm(eHive.BaseRunnable):
    """Normalize a VCF file using bcftools norm"""

    def run(self):
        filepath = self.param_required('filepath')

        self.warning('Analysing file: %s'% filepath)

        file = os.path.split(filepath)[1]
        work_dir = None
        if self.param_is_defined('work_dir'):
            if not os.path.isdir(self.param('work_dir')):
                os.makedirs(self.param('work_dir'))
            work_dir = self.param('work_dir')
        else:
            work_dir = os.path.split(filepath)[0]

        outprefix = "{0}/{1}".format(self.param('work_dir'), file)

        vcfNorm = VcfNormalize(vcf=filepath, bcftools_folder=self.param_required('bcftools_folder'))

        vcf_file = vcfNorm.run_bcftoolsnorm(outprefix=outprefix,
                                            reference=self.param_required('reference'),
                                            multiallelics=self.param('multiallelics'),
                                            type=self.param('type'))

        self.param('out_vcf', vcf_file)

    def write_output(self):
        self.warning('Work is done!')

        self.dataflow({'out_vcf': self.param('out_vcf')}, 1)
