import eHive
import os
from VCFIntegration.SNPTools import SNPTools

class SNPTools_prob2vcf(eHive.BaseRunnable):
    """Run SNPTools prob2vcf on a VCF containing biallelic SNPs"""
    
    def run(self):
        self.warning("Analysing file {0}".format(self.param_required('vcf_file')))

        vcf_i=SNPTools(vcf=self.param_required('vcf_file'),snptools_folder=self.param_required('snptools_folder'))
        chro=self.param_required('chr').rstrip('\n')

        outprefix=os.path.split(self.param_required('outprefix'))[1]
        
        vcf_f=""
        if self.param_is_defined('verbose'):
            vcf_f=vcf_i.run_prob2vcf(probf=self.param_required('probf'),outprefix=outprefix+".{0}".format(chro),
                                              chro=self.param_required('chr'),outdir=self.param_required('work_dir'), verbose=True)
        else:
            vcf_f=vcf_i.run_prob2vcf(probf=self.param_required('probf'),outprefix=outprefix+".{0}".format(chro),
                                              chro=self.param_required('chr'),outdir=self.param_required('work_dir'), verbose=False)
       
        self.param('vcf_f', vcf_f)

    def write_output(self):
        self.warning('Work is done!')

        self.dataflow( {'vcf_f' : self.param('vcf_f') }, 1)


