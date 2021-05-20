import eHive
import os
import re
from VCF.VCFIntegration.Shapeit import Shapeit

class run_Shapeit_convert2vcf(eHive.BaseRunnable):
    """Run SHAPEIT's -convert to generate a VCF file"""

    def fetch_input(self):

        hap_gz = re.sub('\.gz$', '', self.param_required('hap_gz'))
        hap_sample = re.sub('\.sample$', '', self.param_required('hap_sample'))

        if not hap_gz == hap_sample:
            raise Exception("Correct the input prefixes. They are not the same")

        self.param('input_prefix', hap_gz)

    def run(self):

        shapeit_o = Shapeit(shapeit_folder=self.param_required('shapeit_folder'))

        verbose = None
        if self.param_is_defined('verbose'):
            verbose = True
        else:
            verbose = False

        if not os.path.isdir(self.param_required('work_dir')):
            os.makedirs(self.param_required('work_dir'))

        compress = None
        if self.param_is_defined('compress'):
            compress = True
        else:
            compress = False

        outprefix = os.path.split(self.param_required('outprefix'))[1]
        outprefix = "{0}/{1}".format(self.param_required('work_dir'), outprefix)

        outdict = None
        out_vcf = shapeit_o.convert2vcf(input_prefix=self.param('input_prefix'),
                                        output_prefix=outprefix,
                                        compress=compress, verbose=verbose,
                                        logfile='{0}.log'.format(outprefix))

        self.param('out_vcf', out_vcf)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow({'out_vcf': self.param('out_vcf')}, 1)
