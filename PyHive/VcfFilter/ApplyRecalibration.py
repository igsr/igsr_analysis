import eHive
import os
import datetime
from VCFfilter.GATK import GATK

class ApplyRecalibration(eHive.BaseRunnable):
    """run GATK ApplyRecalibration, which is part of the VQSR filtering procedure"""

    def run(self):

        self.warning('Analysing file: %s'% self.param_required('filepath'))

        VcfFilterO = GATK(vcf=self.param_required('filepath'),caller=self.param_required('caller'),gatk_folder=self.param_required('gatk_folder'), reference=self.param_required('reference'), bgzip_folder=self.param('bgzip_folder'), tabix_folder=self.param('tabix_folder'))

        ts_filter_level=None
        if self.param_is_defined('ts_filter_level'):
            ts_filter_level=self.param('ts_filter_level')

        outfile=VcfFilterO.run_applyrecalibration(mode=self.param_required('mode'), recal_file=self.param_required('recal_file'), 
                                                  ts_filter_level=ts_filter_level,tranches_file=self.param_required('tranches_file'), 
                                                  outprefix=self.param_required('filepath'))

        self.param('vcf_filt', outfile)
        self.param('vcf_filt_ix', outfile+".tbi")

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 
            'vcf_filt' : self.param('vcf_filt'),
            'vcf_filt_ix' : self.param('vcf_filt_ix')
        }, 1)

