import eHive
import os

from VcfUtils import VcfUtils

class convertPL2GL(eHive.BaseRunnable):
    """Convert PL fields in the VCF to GL"""

    def run(self):
        filepath=self.param_required('filepath')

        self.warning('Analysing file: %s'% filepath)

        if not os.path.isdir(self.param_required('work_dir')):
            os.makedirs(self.param_required('work_dir'))

        outprefix=os.path.split(self.param_required('outprefix'))[1]

        vcf_object=VcfUtils(vcf=filepath,
                            bcftools_folder=self.param_required('bcftools_folder'))

        vcf_file=vcf_object.convert_PL2GL(outfile=outprefix+".GL.vcf.gz", verbose=True)

        self.param('vcf_file', vcf_file)

    def write_output(self):
        self.warning('Work is done!')

        self.dataflow( {'vcf_file' : self.param('vcf_file') }, 1)
