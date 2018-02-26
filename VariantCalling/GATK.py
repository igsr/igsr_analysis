'''
Created on 22 Feb 2018

@author: ernesto
'''

import os
import subprocess
import pdb

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

    def run_ug(self, outprefix, glm='SNP', compress=True, nt=1, intervals=None, verbose=None, **kwargs):
        
        '''
        Run GATK UnifiedGenotyper

        Parameters
        ----------
        outprefix : str, Required
                    Prefix for output VCF file. i.e. /path/to/file/test
        glm : str, Required
              Genotype likelihoods calculation model to employ -- SNP is the default option, 
              while INDEL is also available for calling indels and BOTH is available for
              calling both together
        compress : boolean, Default= True
                   Compress the output VCF
        nt : int, Optional
             Number of data threads to allocate to UG
        intervals : str, Optional
                    Path to file with genomic intervals to operate with. Also coordinates
                    can be set directly on the command line. For example: chr1:100-200
        verbose : bool, optional
                  if true, then print the command line used for running this program

        Returns
        -------
        A VCF file

        '''

        command = "java -jar "
        if self.gatk_folder:
            command += self.gatk_folder + "/"

        command += "GenomeAnalysisTK.jar -T UnifiedGenotyper -R {0} -I {1} " \
                   "-glm {2} ".format(self.reference, self.bam, glm)

        if intervals:
            command += " -L {0}".format(intervals)

        for k,v in kwargs.items():
            command += " --{0} {1}".format(k,v)

        if compress is True:
            outprefix += ".vcf.gz"
            command += "| {0}/bgzip -c > {1}".format(self.bgzip_folder,outprefix)

        else:
            outprefix += ".vcf"
            command += " -o {0}".format(outprefix)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exp:
            print("Something went wrong while running GATK UG")
            raise Exception(exp.output)

        return outprefix
