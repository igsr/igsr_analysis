import eHive
import os
import pdb
import time
from VCF.VCFfilter.GATK import GATK

class VariantRecalibrator(eHive.BaseRunnable):
    """run GATK VariantRecalibrator, which is part of the VQSR filtering procedure"""

    def run(self):

        self.warning('Analysing file: %s'% self.param_required('filepath'))
        
        vcf = GATK(vcf=self.param_required('filepath'),caller=self.param_required('caller'),gatk_folder=self.param_required('gatk_folder'), reference=self.param_required('reference'))

        optional_params={}
        if self.param_is_defined('annotations'):
            optional_params['annotations']=self.param('annotations')
        if self.param_is_defined('intervals'):
            optional_params['intervals']=self.param('intervals')
        if self.param_is_defined('max_gaussians'):
            optional_params['max_gaussians']=self.param('max_gaussians')
        if self.param_is_defined('tranches'):
            optional_params['tranches']=self.param('tranches')
        
        verbose=False
        if self.param_is_defined('verbose'):
            verbose=True

        log_file=None
        if self.param_is_defined('log_file'):
            log_file="{0}_{1}.log".format(self.param('log_file'),time.strftime("%Y%m%d_%H%M%S"))

        d_out=vcf.run_variantrecalibrator(self.param_required('resources'),mode=self.param_required('mode'), outprefix=self.param_required('filepath'), verbose=verbose, 
                                          log_file=log_file, **optional_params)
        
        self.param('recal_f', d_out['recal_f'])
        self.param('tranches_f', d_out['tranches_f'])

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 
            'recal_f' : self.param('recal_f'),
            'tranches_f' : self.param('tranches_f')
        }, 1)


