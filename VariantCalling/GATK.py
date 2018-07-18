'''
Created on 22 Feb 2018

@author: ernesto
'''

import os
import subprocess
from Utils.RunProgram import RunProgram
from collections import namedtuple
import pdb
import re

class GATK(object):
    '''
    Class to run the different GATK variant callers
    '''


    def __init__(self, bam, reference, gatk_folder, bgzip_folder=None):
        '''
        Constructor

        Class variables
        ---------------
        bam : str, Required
             Path to BAM file used for the variant calling process
        reference : str, Required
             Path to fasta file containing the reference
        gatk_folder : str, Required
                       Path to folder containing GATK's jar file
        bgzip_folder : str, Optional
                       Path to folder containing the bgzip binary
        '''

        if os.path.isfile(bam) is False:
            raise Exception("BAM file does not exist")

        self.bam = bam
        self.reference = reference
        self.gatk_folder = gatk_folder
        self.bgzip_folder = bgzip_folder

    def run_ug(self, outprefix, glm='SNP', compress=True, nt=1, verbose=None, intervals=None, log_file=None, **kwargs):
        
        '''
        Run GATK UnifiedGenotyper

        Parameters
        ----------
        outprefix : str, Required
                    Prefix for output VCF file. i.e. /path/to/file/test
        glm : str, Required
              Genotype likelihoods calculation model to employ -- SNP is the default option, 
              while INDEL is also available for calling indels and BOTH is available for
              calling both together. Default= SNP
        compress : boolean, Default= True
                   Compress the output VCF
        nt : int, Optional
             Number of data threads to allocate to UG
        intervals : : list, Optional
                    List in which each of the elements is a path to file with genomic intervals to operate with. Also coordinates
                    can be set directly on the command line. For example: ['chr1:100-200', 'chr2:200-300']. If the list contains
                    more than one interval, then it is useful to set the --interval_set_rule option
        verbose : bool, optional
                  if true, then print the command line used for running this program
        alleles: str, Optional
                 Path to VCF.
                 When --genotyping_mode is set to
                 GENOTYPE_GIVEN_ALLELES mode, the caller will genotype the samples
                 using only the alleles provide in this callset
        genotyping_mode: str, Optional
                         Specifies how to determine the alternate alleles to use for genotyping
                         Possible values are: DISCOVERY, GENOTYPE_GIVEN_ALLELES
        output_mode: str, Optional
                     Which type of calls we should output.
                     Possible values are: EMIT_VARIANTS_ONLY, EMIT_ALL_CONFIDENT_SITES, EMIT_ALL_SITES
                     Default: EMIT_VARIANTS_ONLY
        log_file : str, Optional
                   Path to file that will used for logging the GATK stderr and stdout

        Returns
        -------
        A VCF file

        '''

        Arg = namedtuple('Argument', 'option value')

        arguments=[Arg('-T','UnifiedGenotyper'), Arg('-R',self.reference), Arg('-I',self.bam), Arg('-glm',glm), Arg('-nt',nt)]

        if intervals is not None:
            for i in intervals:
                arguments.append(Arg('--intervals',i))
        
        for k,v in kwargs.items():
            if v is not None: arguments.append(Arg(" --{0}".format(k),v))

        pipelist=None
        if compress is True:
            outprefix += ".vcf.gz"
            compressRunner=RunProgram(path=self.bgzip_folder,program='bgzip',parameters=[ '-c', '>', outprefix])
            pipelist=[compressRunner]
        else:
            outprefix += ".vcf"
            arguments.append(Arg('-o',outprefix))

        runner=RunProgram(program='java -jar {0}/GenomeAnalysisTK.jar'.format(self.gatk_folder), args=arguments, downpipe=pipelist, log_file=log_file)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))
            
        stdout,stderr=runner.run_popen()

        return outprefix

    def run_hc(self, outprefix, compress=True, verbose=None, log_file=None, intervals=None, **kwargs):
        '''
        Run GATK HaplotypeCaller

        Parameters
        ----------
        outprefix : str, Required
                    Prefix for output VCF file. i.e. /path/to/file/test
        compress : boolean, Default= True
                   Compress the output VCF
        num_cpu_threads_per_data_thread : int, Optional
                   controls the number of CPU threads allocated to each data thread
        intervals : list, Optional
                    List in which each of the elements is a path to file with genomic intervals to operate with. Also coordinates
                    can be set directly on the command line. For example: ['chr1:100-200', 'chr2:200-300']. If the list contains 
                    more than one interval, then it is useful to set the --interval_set_rule option
        standard_min_confidence_threshold_for_calling : int, Optional
                                                        The minimum phred-scaled confidence threshold at which variants should be called
                                                        Default: 10
        genotyping_mode: str, Optional
                         Specifies how to determine the alternate alleles to use for genotyping
                         Possible values are: DISCOVERY, GENOTYPE_GIVEN_ALLELES
        alleles: str, Optional
                 Path to VCF.
                 When --genotyping_mode is set to
                 GENOTYPE_GIVEN_ALLELES mode, the caller will genotype the samples
                 using only the alleles provide in this callset
        emitRefConfidence: str, Optional
                           Mode for emitting reference confidence scores
                           Possible values are: NONE, BP_RESOLUTION, GVCF
        verbose : bool, optional
                  if true, then print the command line used for running this program
        log_file : str, Optional
                   Path to file that will used for logging the GATK stderr and stdout

        Returns
        -------
        A VCF file

        '''
        Arg = namedtuple('Argument', 'option value')

        arguments=[Arg('-T','HaplotypeCaller'), Arg('-R',self.reference), Arg('-I',self.bam)]
        
        if intervals is not None:
            for i in intervals:
                arguments.append(Arg('--intervals',i)) 

        for k,v in kwargs.items():
            if v is not None: arguments.append(Arg("--{0}".format(k),v))

        pipelist=None
        if compress is True:
            outprefix += ".vcf.gz"
            compressRunner=RunProgram(path=self.bgzip_folder,program='bgzip',parameters=[ '-c', '>', outprefix])
            pipelist=[compressRunner]
        else:
            outprefix += ".vcf"
            arguments.append(Arg('-o',outprefix))
            
        runner=RunProgram(program='java -jar {0}/GenomeAnalysisTK.jar'.format(self.gatk_folder), downpipe=pipelist, args=arguments, log_file=log_file)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout,stderr,is_error=runner.run_popen(raise_exc=False)

        if is_error is True:
            '''
            This piece of code is necessary as GATK crashes when the intersection between the genomic chunk and the alleles passed in the VCF 
            are calculated and there are no sites.

            If that's the case then GATK  will be run without the interval intersection
            '''
            patt = re.compile('##### ERROR MESSAGE: Bad input: The INTERSECTION of your -L options produced no intervals.')
            lines=stderr.split('\n')
            interval_error_seen=False
            for l in lines:
                m = patt.match(l)
                if m:
                    interval_error_seen=True
                    alleles=([arg.value for arg in arguments if arg.option == '--alleles'])[0]
                    for k,i in enumerate(arguments):
                        if i.option == '--intervals' and i.value == alleles:
                            del arguments[k]
                        elif i.option == '--interval_set_rule':
                            del arguments[k]                    
            if interval_error_seen is False:
                 raise Exception(stderr)
            elif interval_error_seen is True:
                runner=RunProgram(program='java -jar {0}/GenomeAnalysisTK.jar'.format(self.gatk_folder), downpipe=pipelist, args=arguments, log_file=log_file)
                if verbose is True:
                    print("Command line is: {0}".format(runner.cmd_line))
                      
                stdout,stderr,is_error=runner.run_popen(raise_exc=True)

        return outprefix
