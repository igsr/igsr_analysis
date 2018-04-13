import eHive
import subprocess
import os
import pdb
import sys
import glob

from GATK import GATK

class GATK_UG(eHive.BaseRunnable):
    '''
    Run GATK UnifiedGenotyper on a BAM file/s

    '''
    
    def run(self):

        if not os.path.isdir(self.param_required('work_dir')):
            os.makedirs(self.param_required('work_dir'))

        outprefix=os.path.split(self.param_required('outprefix'))[1]

        chrom=self.param_required('chunk')[0]
        start=None
        # increment by 1 if start=0, as GATK does not accept coords <1
        if self.param_required('chunk')[1]==0:
            start=1
        else:
            start=self.param_required('chunk')[1]
        end=self.param_required('chunk')[2]
        
        outfile="{0}/{1}.ug.{2}_{3}_{4}".format(self.param_required('work_dir'), 
                                                outprefix, 
                                                chrom,
                                                start,
                                                end)

        gatk_object=GATK(bam=self.param_required('bamlist'), 
                         reference=self.param_required('reference'), 
                         gatk_folder=self.param_required('gatk_folder'),
                         bgzip_folder=self.param_required('bgzip_folder'))

        intervals=None
        
        if self.param_is_defined('chunk'):
            intervals="{0}:{1}-{2}".format(chrom, 
                                           start,
                                           end)

        alleles=None
        if self.param_is_defined('alleles'):
            alleles=self.param('alleles')

        genotyping_mode=None
        if self.param_is_defined('genotyping_mode'):
            genotyping_mode=self.param('genotyping_mode')

        nt=None
        if self.param_is_defined('threads'):
            nt=self.param('threads')

        max_deletion_fraction=0.05
        if self.param_is_defined('max_deletion_fraction'):
            max_deletion_fraction=self.param('max_deletion_fraction')

        dcov=250
        if self.param_is_defined('dcov'):
            dcov=self.param('dcov')

        outfile=gatk_object.run_ug(outprefix=outfile,
                                   glm=self.param_required('glm'),
                                   output_mode=self.param_required('output_mode'),
                                   downsample_to_coverage=dcov,
                                   alleles=alleles, genotyping_mode=genotyping_mode,
                                   intervals=intervals, nt=nt, 
                                   max_deletion_fraction=max_deletion_fraction)

        self.param('out_vcf', outfile)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'out_vcf' : self.param('out_vcf') }, 1)




