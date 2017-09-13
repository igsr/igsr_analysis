import eHive
import os
import datetime
from VcfFilter import VcfFilter

class ApplyRecalibration(eHive.BaseRunnable):
    """run GATK ApplyRecalibration, which is part of the VQSR filtering procedure"""

    def run(self):

        self.warning('Analysing file: %s'% self.param_required('filepath'))

        VcfFilterO = VcfFilter(vcf=self.param_required('filepath'),caller=self.param_required('caller'),gatk_folder=self.param_required('gatk_folder'), reference=self.param_required('reference'), bgzip_folder=self.param('bgzip_folder'), tabix_folder=self.param('tabix_folder'))

        outfile=VcfFilterO.run_applyrecalibration(mode=self.param_required('mode'), recal_file=self.param_required('recal_file'), tranches_file=self.param_required('tranches_file'), outprefix=self.param_required('filepath'))

        self.param('vcf_filt', outfile)
        self.param('vcf_filt_ix', outfile+".tbi")

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'vcf_filt' : self.param('vcf_filt') }, 1)
        self.dataflow( { 'vcf_filt_ix' : self.param('vcf_filt_ix') }, 1)
