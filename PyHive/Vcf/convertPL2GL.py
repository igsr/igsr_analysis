import eHive
import os

from VCF.VcfUtils import VcfUtils

class convertPL2GL(eHive.BaseRunnable):
    """Convert PL fields in the VCF to GL"""

    def run(self):
        filepath = self.param_required('filepath')

        self.warning('Analysing file: %s'% filepath)

        if not os.path.isdir(self.param_required('work_dir')):
            os.makedirs(self.param_required('work_dir'))

        threads = 1
        if self.param_is_defined('threads'):
            threads = self.param('threads')

        outprefix = os.path.split(self.param_required('outprefix'))[1]

        outfile = "{0}/{1}.GL.vcf.gz".format(self.param_required('work_dir'), outprefix)

        vcf_object = VcfUtils(vcf=filepath,
                              bcftools_folder=self.param_required('bcftools_folder'))

        vcf_file = vcf_object.convert_PL2GL(outfile, threads=threads, verbose=True)

        self.param('out_vcf', vcf_file)

    def write_output(self):
        self.warning('Work is done!')

        self.dataflow({'out_vcf': self.param('out_vcf')}, 1)
