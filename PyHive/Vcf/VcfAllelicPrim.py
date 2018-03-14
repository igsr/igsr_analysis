import eHive
import os
import datetime
from VcfNormalize import VcfNormalize

class VcfAllelicPrim(eHive.BaseRunnable):
    """Run vcfallelicprimitives on a VCF file"""

    def run(self):
        filepath=self.param_required('filepath')

        self.warning('Analysing file: %s'% filepath)

        file=os.path.split(filepath)[1]
        work_dir=None
        if self.param_is_defined('work_dir'):
            work_dir=self.param_required('work_dir')
        else:
            work_dir=os.path.split(filepath)[0]

        vcfNorm = VcfNormalize(vcf=filepath, vcflib_folder=self.param_required('vcflib_folder'), bgzip_folder=self.param_required('bgzip_folder'))

        downstream_pipe=None

        vcf_file=""
        if self.param_is_defined('compress'):
            vcf_file=vcfNorm.run_vcfallelicprimitives(outprefix=file,compress=True,downstream_pipe=downstream_pipe,outdir=work_dir)
        else:
            vcf_file=vcfNorm.run_vcfallelicprimitives(outprefix=file,compress=False,downstream_pipe=downstream_pipe,outdir=work_dir)

        self.param('out_vcf', vcf_file)

    def write_output(self):
        self.warning('Work is done!')

        self.dataflow( {'out_vcf' : self.param('out_vcf') }, 1)
