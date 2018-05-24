import eHive
import subprocess
import os
import pdb
import sys
import glob

from VariantCalling import GATK

class GATK_UG(eHive.BaseRunnable):
    '''
    Run GATK UnifiedGenotyper on a BAM file/s

    eHive runnable for GATK UG caller

    Parameters
    ----------
    work_dir : str, Required
         Path to the working directory
    outprefix : str, Required
         String used as prefix for the output file
    chunk: str, Required
         Interval used by GATK UG for the variant calling
    bamlist: str, Required
         Path to BAM file used for the variant calling
    reference: str, Required
         Path to the Fasta file used to generate the BAM alignment file
    gatk_folder: str, Required
         Path to folder containing the GATK jar file
    bgzip_folder: str, Required
         Path to folder containing the Bgzip binary
    glm: str, Required
         Genotype likelihoods calculation model to employ.
         Possible values are: SNP, INDEL, BOTH
    output_mode: str, Required
         Which type of calls we should output.
         Possible values are: EMIT_VARIANTS_ONLY, EMIT_ALL_CONFIDENT_SITES, EMIT_ALL_SITES
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
    max_deletion_fraction: float, Optional
         Maximum fraction of reads with deletions spanning this locus for it to be callable
         Default value=0.05. Note: To disable the use of this parameter, set its value to >1.
    dcov: int, Optional
         Target coverage threshold for downsampling to coverage.
         Default value=250

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
        
        outfile="{0}/{1}.ug.{2}_{3}_{4}".format(self.param_required('work_dir'), 
                                                outprefix, 
                                                chrom,
                                                start,
                                                end)

        gatk_object=GATK(bam=self.param_required('bamlist'), 
                         reference=self.param_required('reference'), 
                         gatk_folder=self.param_required('gatk_folder'),
                         bgzip_folder=self.param_required('bgzip_folder'))

        intervals="{0}:{1}-{2}".format(chrom, 
                                       start,
                                       end)

        alleles=None
        if self.param_is_defined('alleles'):
            alleles=self.param('alleles')

        genotyping_mode='DISCOVERY'
        if self.param_is_defined('genotyping_mode'):
            genotyping_mode=self.param('genotyping_mode')

        nt=1
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
                                   alleles=alleles, 
                                   genotyping_mode=genotyping_mode,
                                   intervals=intervals, nt=nt, 
                                   max_deletion_fraction=max_deletion_fraction)

        self.param('out_vcf', outfile)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'out_vcf' : self.param('out_vcf') }, 1)




