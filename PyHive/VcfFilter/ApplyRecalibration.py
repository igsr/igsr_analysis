import eHive
import os
import time
from VCF.VCFfilter.GATK import GATK

class ApplyRecalibration(eHive.BaseRunnable):
    """run GATK ApplyRecalibration, which is part of the VQSR filtering procedure"""

    def run(self):

        self.warning('Analysing file: %s'% self.param_required('filepath'))

        VcfFilterO = GATK(vcf=self.param_required('filepath'),
                          caller=self.param_required('caller'),
                          gatk_folder=self.param_required('gatk_folder'),
                          reference=self.param_required('reference'),
                          bgzip_folder=self.param('bgzip_folder'),
                          tabix_folder=self.param('tabix_folder'),
                          tmp_dir=self.param('tmp_dir'))

        ts_filter_level = None
        if self.param_is_defined('ts_filter_level'):
            ts_filter_level = self.param('ts_filter_level')

        outprefix = None
        filename = os.path.split(self.param_required('filepath'))[1]
        work_dir = None
        if self.param_is_defined('work_dir'):
            work_dir = self.param_required('work_dir')
        else:
            work_dir = os.path.split(self.param_required('filepath'))[0]
        outprefix = "{0}/{1}".format(work_dir, filename)

        threads = 1
        if self.param_is_defined('threads'):
            threads = self.param('threads')

        log_file = None
        if self.param_is_defined('log_file'):
            log_file = "{0}_{1}.log".format(self.param('log_file'), time.strftime("%Y%m%d_%H%M%S"))

        outfile = VcfFilterO.run_applyrecalibration(mode=self.param_required('mode'),
                                                    recal_file=self.param_required('recal_file'),
                                                    ts_filter_level=ts_filter_level,
                                                    tranches_file=self.param_required('tranches_file'),
                                                    num_threads=threads, outprefix=outprefix,
                                                    log_file=log_file)

        self.param('vcf_filt', outfile)
        self.param('vcf_filt_ix', outfile+".tbi")

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow({
            'vcf_filt': self.param('vcf_filt'),
            'vcf_filt_ix': self.param('vcf_filt_ix')
        }, 1)
