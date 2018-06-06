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

    def run_ug(self, outprefix, glm='SNP', compress=True, nt=1, verbose=None, **kwargs):
        
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
        intervals : str, Optional
                    Path to file with genomic intervals to operate with. Also coordinates
                    can be set directly on the command line. For example: chr1:100-200
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

        Returns
        -------
        A VCF file

        '''

        Arg = namedtuple('Argument', 'option value')

        arguments=[Arg('-T','UnifiedGenotyper'), Arg('-R',self.reference), Arg('-I',self.bam), Arg('-glm',glm), Arg('-nt',nt)]

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

        runner=RunProgram(program='java -jar {0}/GenomeAnalysisTK.jar'.format(self.gatk_folder), args=arguments, downpipe=pipelist)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))
            
        stdout,stderr=runner.run_popen()

        return outprefix
