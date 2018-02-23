import eHive
import os
import datetime
from VcfNormalize import VcfNormalize

class NormVcf(eHive.BaseRunnable):
    """Normalize a VCF file using vt normalize"""

    def param_defaults(self):
        return {
        }

    def run(self):
        filepath=self.param_required('filepath')

        self.warning('Analysing file: %s'% filepath)

        file=os.path.split(filepath)[1]
        work_dir=None
        if self.param_is_defined('work_dir'):
            work_dir=self.param_required('work_dir')
        else:
            work_dir=os.path.split(filepath)[0]

        vcfNorm = VcfNormalize(vcf=filepath, vt_folder=self.param_required('vt_folder'), bgzip_folder=self.param_required('bgzip_folder'))
        vcf_file=""
        if self.param_is_defined('compress'):
            vcf_file=vcfNorm.run_vtnormalize(outprefix=file,reference=self.param_required('reference'),compress=True,outdir=work_dir)
        else:
            vcf_file=vcfNorm.run_vtnormalize(outprefix=file,reference=self.param_required('reference'),compress=False,outdir=work_dir)

        self.param('vcf_file', vcf_file)

    def write_output(self):
        self.warning('Work is done!')

        self.dataflow( {'vcf_file' : self.param('vcf_file') }, 1)
