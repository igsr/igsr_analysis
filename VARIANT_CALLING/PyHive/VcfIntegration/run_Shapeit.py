import eHive
import os
import pdb
from VCFIntegration.Shapeit import Shapeit

class run_Shapeit(eHive.BaseRunnable):
    """Run SHAPEIT"""
    
    def fetch_input(self):
        if self.param_is_defined('chunk'):
            self.param('inputfrom', self.param('chunk')[1])
            self.param('inputto', self.param('chunk')[2])

    def run(self):
    
        vcf_g=Shapeit(shapeit_folder=self.param_required('shapeit_folder'))

        verbose=None
        if self.param_is_defined('verbose'):
            verbose=True
        else:
            verbose=False

        outprefix=os.path.split(self.param_required('outprefix'))[1]
        outprefix="{0}/{1}".format(self.param_required('work_dir'),outprefix)
            
        options_dict={}
        
        if self.param_is_defined('inputthr'):
            options_dict['input-thr']=self.param('inputthr')
        if self.param_is_defined('thread'):
            options_dict['thread']=self.param('thread')
        if self.param_is_defined('window'):
            options_dict['window']=self.param('window')
        if self.param_is_defined('statesrandom'):
            options_dict['states-random']=self.param('statesrandom')
        if self.param_is_defined('burn'):
            options_dict['burn']=self.param('burn')
        if self.param_is_defined('run'):
            options_dict['run']=self.param('run')
        if self.param_is_defined('prune'):
            options_dict['prune']=self.param('prune')
        if self.param_is_defined('main'):
            options_dict['main']=self.param('main')
        if self.param_is_defined('inputfrom'):
            outprefix += ".{0}".format(self.param('inputfrom'))
            options_dict['input-from']=self.param('inputfrom')
        if self.param_is_defined('inputto'):
            outprefix += ".{0}".format(self.param('inputto'))
            options_dict['input-to']=self.param('inputto')

        duohmm=False
        if self.param_is_defined('duohmm'):
            duohmm=True

        shapeit_o=Shapeit(shapeit_folder = self.param_required('shapeit_folder'))

        outdict=None
        if self.param_is_defined('input_gen'):       
            outdict=shapeit_o.run_shapeit(input_gen= self.param_required('input_gen'),
                                          input_init= self.param_required('input_init'),
                                          input_scaffold= self.param_required('input_scaffold'),
                                          output_prefix= outprefix,
                                          duohmm= duohmm,
                                          verbose=verbose,
                                          **options_dict)
        elif self.param_is_defined('input_bed'):
            outdict=shapeit_o.run_shapeit(input_bed= self.param_required('input_bed'),
                                          duohmm= duohmm,
                                          output_prefix= outprefix,
                                          verbose=verbose,
                                          **options_dict)
        self.param('outdict', outdict)
       
    def write_output(self):
        self.warning('Work is done!')

        outdict= self.param('outdict')
        self.dataflow( {'hap_gz' : outdict['hap_gz'] }, 1)
        self.dataflow( {'hap_sample' : outdict['hap_sample'] }, 1)




