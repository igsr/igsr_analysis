import eHive
import os
from VcfGenotype import VcfGenotype

class SNPTools_prob2vcf(eHive.BaseRunnable):
    """Run SNPTools prob2vcf on a VCF containing biallelic SNPs"""
    
    def run(self):
        self.warning("Analysing file {0}".format(self.param_required('vcf_file')))

        vcf_g=VcfGenotype(vcf=self.param_required('vcf_file'),snptools_folder=self.param_required('snptools_folder'))
        chro=self.param_required('chr').rstrip('\n')

        outprefix=self.param_required('outprefix')
        
        vcf_f=""
        if self.param_is_defined('verbose'):
            vcf_f=vcf_g.run_snptools_prob2vcf(probf=self.param_required('probf'),outprefix=outprefix+"_{0}".format(chro),
                                              chro=self.param_required('chr'),outdir=self.param_required('work_dir'), verbose=True)
        else:
            vcf_f=vcf_g.run_snptools_prob2vcf(probf=self.param_required('probf'),outprefix=outprefix+"_{0}".format(chro),
                                              chro=self.param_required('chr'),outdir=self.param_required('work_dir'), verbose=False)
       
        self.param('vcf_f', vcf_f)

    def write_output(self):
        self.warning('Work is done!')

        self.dataflow( {'vcf_f' : self.param('vcf_f') }, 1)


