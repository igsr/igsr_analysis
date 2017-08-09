import eHive
from VcfGenotype import VcfGenotype

class SNPTools_poprob(eHive.BaseRunnable):
    """Run SNPTools poprob on a VCF containing biallelic SNPs"""
    
    def run(self):
       vcf_g=VcfGenotype(vcf=self.param_required('vcf_file'),snptools_folder=self.param_required('snptools_folder'))
       
       prob_f=vcf_g.run_snptools_poprob(outprefix=self.param_required('outprefix'), rawlist=self.param_required('rawlist'), outdir=self.param_required('work_dir'))
       
       self.param('prob_f', prob_f)

    def write_output(self):
        self.warning('Work is done!')

        self.dataflow( {'prob_f' : self.param('prob_f') }, 1)


