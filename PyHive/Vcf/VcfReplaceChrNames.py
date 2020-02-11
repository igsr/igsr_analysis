import eHive
import os

from VCF.VcfUtils import VcfUtils

class VcfReplaceChrNames(eHive.BaseRunnable):
    """Modify the style (UCSC or Ensembl) of the chr names in a VCF file"""

    def run(self):
        filepath = self.param_required('filepath')

        self.warning('Analysing file: %s'% filepath)

        basename = os.path.split(filepath)[1]
        work_dir = None
        if self.param_is_defined('work_dir'):
            if not os.path.isdir(self.param('work_dir')):
                os.makedirs(self.param('work_dir'))
            work_dir = self.param('work_dir')
        else:
            work_dir = os.path.split(filepath)[0]

        vcf_object = VcfUtils(vcf=filepath, bgzip_folder=self.param('bgzip_folder'))

        outprefix = "{0}/{1}".format(work_dir, basename+".{0}".
                                     format(self.param_required('chr_types')))

        outfile = vcf_object.rename_chros(chr_types=self.param_required('chr_types'),
                                          outfile=outprefix+".vcf.gz")

        self.param('vcf_f', outfile)

    def write_output(self):
        self.warning('Work is done!')

        self.dataflow({'vcf_f': self.param('vcf_f')}, 1)
