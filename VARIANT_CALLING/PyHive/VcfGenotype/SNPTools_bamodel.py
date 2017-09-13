import eHive
from VcfGenotype import VcfGenotype

class SNPTools_bamodel(eHive.BaseRunnable):
    """Run SNPTools bamodel on a VCF containing biallelic SNPs"""
    
    def run(self):
       vcf_g=VcfGenotype(vcf=self.param_required('vcf_file'),snptools_folder=self.param_required('snptools_folder'))

       raw_f=""
       if self.param_is_defined('verbose'):
           raw_f=vcf_g.run_snptools_bamodel(sample=self.param_required('sample'), bamfiles=self.param_required('bamlist'), outdir=self.param_required('work_dir'), verbose=True)
       else:
           raw_f=vcf_g.run_snptools_bamodel(sample=self.param_required('sample'), bamfiles=self.param_required('bamlist'), outdir=self.param_required('work_dir'), verbose=False)

       self.param('raw_f', raw_f)
       
    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( {'raw_f' : self.param('raw_f') }, 1)



