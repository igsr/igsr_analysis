import eHive
import os
from VCF.VcfNormalize import VcfNormalize

class VariantsToAllelicPrimitives(eHive.BaseRunnable):
    """Run GATK VariantsToAllelicPrimitives on a VCF file"""

    def run(self):
        filepath = self.param_required('filepath')

        self.warning('Analysing file: %s'% filepath)

        file = os.path.split(filepath)[1]
        work_dir = None
        if self.param_is_defined('work_dir'):
            work_dir = self.param_required('work_dir')
        else:
            work_dir = os.path.split(filepath)[0]

        vcfNorm = VcfNormalize(vcf=filepath,
                               gatk_folder=self.param_required('gatk_folder'),
                               bgzip_folder=self.param_required('bgzip_folder'))
        vcf_file = ""

        if self.param_is_defined('compress'):
            compress = None
            if self.param('compress') == 'True':
                compress = True
            elif self.param('compress') == 'False':
                compress = False
            else:
                raise Exception('compress parameter is not valid: {0}'.
                                format(self.param('compress')))
            vcf_file = vcfNorm.run_gatk_VariantsToAllelicPrimitives(outprefix=file,
                                                                    reference=self.param_required('reference'),
                                                                    outdir=work_dir,
                                                                    compress=compress)
        else:
            vcf_file = vcfNorm.run_gatk_VariantsToAllelicPrimitives(outprefix=file,
                                                                    reference=self.param_required('reference'),
                                                                    outdir=work_dir)

        self.param('vcf_file', vcf_file)

    def write_output(self):
        self.warning('Work is done!')

        self.dataflow({'vcf_file': self.param('vcf_file')}, 1)
