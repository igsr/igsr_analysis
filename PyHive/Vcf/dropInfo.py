import eHive
import os
import pdb

from VcfUtils import VcfUtils

class dropInfo(eHive.BaseRunnable):
    """drop INFO annotation from a VCF file"""

    def run(self):
        filepath=self.param_required('filepath')

        self.warning('Analysing file: %s'% filepath)

        if not os.path.isdir(self.param_required('work_dir')):
            os.makedirs(self.param_required('work_dir'))

        outprefix=os.path.split(self.param_required('outprefix'))[1]

        outfile = "{0}/{1}.noINFO.vcf.gz".format(self.param_required('work_dir'), outprefix)

        vcf_object=VcfUtils(vcf=filepath,
                            bcftools_folder=self.param_required('bcftools_folder'))

        vcf_file=vcf_object.drop_info(outfile, verbose=True)

        self.param('out_vcf', vcf_file)

    def write_output(self):
        self.warning('Work is done!')

        self.dataflow( {'out_vcf' : self.param('out_vcf') }, 1)
