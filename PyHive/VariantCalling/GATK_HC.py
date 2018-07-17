import eHive
import subprocess
import os
import pdb
import sys
import glob
import time

from VariantCalling import GATK

class GATK_HC(eHive.BaseRunnable):
    '''
    Run GATK HaplotypeCaller on a BAM file/s

    eHive runnable for GATK HC caller

    Parameters
    ----------
    work_dir : str, Required
         Path to the working directory
    outprefix : str, Required
         String used as prefix for the output file
    chunk: str, Required
         Interval used by GATK HC for the variant calling
    bamlist: str, Required
         Path to BAM file used for the variant calling
    reference: str, Required
         Path to the Fasta file used to generate the BAM alignment file
    gatk_folder: str, Required
         Path to folder containing the GATK jar file
    bgzip_folder: str, Required
         Path to folder containing the Bgzip binary
    alleles: str, Optional
         Path to VCF. 
         When --genotyping_mode is set to 
         GENOTYPE_GIVEN_ALLELES mode, the caller will genotype the samples 
         using only the alleles provide in this callset
    genotyping_mode: str, Optional
         Specifies how to determine the alternate alleles to use for genotyping
         Possible values are: DISCOVERY, GENOTYPE_GIVEN_ALLELES
    threads: int, Optional
         Number of CPUs used by the caller
         Default=1
    log_file: str, Optional
            Path to log file used to log the GATK HC stderr

    Returns
    -------

    Path to VCF file
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
        
        outfile="{0}/{1}.hc.{2}_{3}_{4}".format(self.param_required('work_dir'), 
                                                outprefix, 
                                                chrom,
                                                start,
                                                end)

        gatk_object=GATK(bam=self.param_required('bamlist'), 
                         reference=self.param_required('reference'), 
                         gatk_folder=self.param_required('gatk_folder'),
                         bgzip_folder=self.param_required('bgzip_folder'))

        intervals=["{0}:{1}-{2}".format(chrom, 
                                       start,
                                       end)]

        interval_set_rule="UNION"
        alleles=None
        if self.param_is_defined('alleles'):
            alleles=self.param('alleles')
            intervals.append(alleles)
            interval_set_rule="INTERSECTION"

        genotyping_mode='DISCOVERY'
        if self.param_is_defined('genotyping_mode'):
            genotyping_mode=self.param('genotyping_mode')

        nt=1
        if self.param_is_defined('threads'):
            nt=self.param('threads')

        standard_min_confidence_threshold_for_calling=10
        if self.param_is_defined('standard_min_confidence_threshold_for_calling'):
            standard_min_confidence_threshold_for_calling=self.param('standard_min_confidence_threshold_for_calling')

        log_file=None
        if self.param_is_defined('log_file'):
            log_file="{0}_{1}.log".format(self.param('log_file'),time.strftime("%Y%m%d_%H%M%S"))

        verbose=None
        if self.param_is_defined('verbose'):
            verbose=True
        else:
            verbose=False

        outfile=gatk_object.run_hc(outprefix=outfile,
                                   alleles=alleles, 
                                   genotyping_mode=genotyping_mode,
                                   intervals=intervals,
                                   interval_set_rule=interval_set_rule,
                                   num_cpu_threads_per_data_thread=nt,
                                   standard_min_confidence_threshold_for_calling=standard_min_confidence_threshold_for_calling,
                                   log_file=log_file, verbose=verbose)

        self.param('out_vcf', outfile)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'out_vcf' : self.param('out_vcf') }, 1)




