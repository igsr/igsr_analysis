import eHive
from VcfGenotype import VcfGenotype

class run_Beagle(eHive.BaseRunnable):
    """Run Beagle on a VCF containing Genotype likelihoods"""
    
    def run(self):
    
        vcf_g=VcfGenotype(vcf=self.param_required('vcf_file'),beagle_folder=self.param_required('beagle_folder'))

        vcf_f=""
        if self.param_is_defined('verbose'):
            vcf_f=vcf_g.run_beagle(outprefix=self.param_required('outprefix'), outdir=self.param_required('work_dir'), window=self.param_required('window'), 
                                   overlap=self.param_required('overlap'), niterations=self.param_required('niterations'), nthreads=self.param_required('nthreads'), verbose=True)
        else:
            vcf_f=vcf_g.run_beagle(outprefix=self.param_required('outprefix'), outdir=self.param_required('work_dir'), window=self.param_required('window'),
                                   overlap=self.param_required('overlap'), niterations=self.param_required('niterations'), nthreads=self.param_required('nthreads'), verbose=False)

        self.param('vcf_f', vcf_f)
       
    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( {'vcf_f' : self.param('vcf_f') }, 1)



