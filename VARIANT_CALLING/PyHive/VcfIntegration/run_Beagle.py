import eHive
import os
from VcfIntegration import VcfIntegration

class run_Beagle(eHive.BaseRunnable):
    """Run Beagle on a VCF containing Genotype likelihoods"""
    
    def run(self):
    
        vcf_g=VcfIntegration(vcf=self.param_required('vcf_file'),beagle_folder=self.param_required('beagle_folder'))

        verbose=None
        if self.param_is_defined('verbose'):
            verbose=True
        else:
            verbose=False

        outprefix=os.path.split(self.param_required('outprefix'))[1]

        options={}
        if self.param_is_defined('window'):
            options['window']=self.param('window')
        if self.param_is_defined('overlap'):
            options['overlap']=self.param('overlap')
        if self.param_is_defined('niterations'):
            options['niterations']=self.param('niterations')
        if self.param_is_defined('nthreads'):
            options['nthreads']=self.param('nthreads')
        
        vcf_f=vcf_g.run_beagle(outprefix=outprefix, outdir=self.param_required('work_dir'), **options, 
                               region=self.param_required('region_chunk'), verbose=verbose)

        self.param('vcf_f', vcf_f)
       
    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( {'vcf_f' : self.param('vcf_f') }, 1)



