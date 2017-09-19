'''
Created on 03 Aug 2017

@author: ernesto
'''

import os
import subprocess
import pdb

class VcfUtils(object):
    '''
    Class to represent a misc of actions that can be done on a single or multiple VCF files
    '''


    def __init__(self, vcf=None, vcflist=None, bcftools_folder=None, gatk_folder=None, java_folder=None):
        '''
        Constructor

        Class variables
        ---------------
        vcf : str, optional
             Path to gzipped vcf file
        vcflist : list, optional
             List of dicts containing setname:vcf_paths (keys:values) pairs
        bcftools_folder : str, Optional
                         Path to folder containing the bcftools binary
        gatk_folder : str, Optional
                      Path to folder containing the jar file
        java_folder : str, Optional
                      Path to folder containing the java binary

        Imp: Either 'vcf' or 'vcflist' variables should be initialized
        '''

        if not vcf and not vcflist:
            raise Exception("Either a vcf file or a list of vcf files should be used\
                             to initialize this class")

        if vcf is not None:
            if os.path.isfile(vcf) is False:
                raise Exception("File does not exist")

        self.vcf = vcf
        self.vcflist = vcflist
        self.bcftools_folder = bcftools_folder
        self.gatk_folder = gatk_folder
        self.java_folder = java_folder

    def reheader(self, newheader, outprefix, samplefile=None, verbose=False):
        '''
        Modifiy the VCF's header with the newheader

        Parameters
        ----------
        newheader : string, required
                   Path to the file containing the new header
        outprefix : string, required
                    prefix for output files
        samplefile : string, optional
                     Path to the file with the sample names that will included
                     in the new header
        verbose : bool, optional
                  increase the verbosity, default=False

        Returns
        -------
        Path to the VCF with the modified header
        '''

        command = ""
        if self.bcftools_folder:
            command += self.bcftools_folder + "/"

        outfile=outprefix+".reheaded.vcf.gz"

        if samplefile is not None:
            command += "bcftools reheader -h {0} -o {1} -s {2} {3}".format(newheader, outfile, samplefile, self.vcf)
        else:
            command += "bcftools reheader -h {0} -o {1} {2}".format(newheader, outfile, self.vcf)

        if verbose is True:
            print("Command is: %s" % command)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            print("Cmd used was {0}".format(exc.cmd))
            raise Exception(exc.output)

        return outfile

    def combine(self, labels, reference, outprefix, compress=False, outdir=None, ginterval=None, 
                genotypemergeoption=None, options=None):
        '''
        Combine VCFs using GATK's CombineVariants into a single VCF

        Parameters
        ----------
        labels : list, required
                 List of labels used for each of the VCFs in self.vcflist. The order of the labels
                 should be the same that the VCFs in the list
        reference : str, required
                    Path to Fasta file with reference
        outprefix : str, required
                    prefix used for output file
        compress : boolean, optional
                   Compress the output VCF with bgzip. Default=False
        outdir : str, optional
                 Path to folder used to write the results to
        ginterval : str, optional
                    Genomic interval used to restrict the analysis. i.e. chr20:1000-2000
        genotypemergeoption : str, optional
                    Determines how we should merge genotype records for samples shared across the ROD files. 
                    Possible values are: UNIQUIFY, PRIORITIZE, UNSORTED, REQUIRE_UNIQUE
        options : list, optional
                   List of options. i.e. ['-env','--filteredAreUncalled']
    
        Returns
        -------
        Path to the merged VCF
        '''

        variants_str=""
        for path, label in zip(self.vcflist, labels):
            if os.path.isfile(path) == False:
                print("Error reading from {0}".format(path))
                raise Exception("File does not exist")
            variants_str+="-V:{0} {1} ".format(label,path)
        
        command = ""
        if self.java_folder:
            command += "{0}/".format(self.java_folder)

        command += "java -jar "

        if self.gatk_folder:
            command += "{0}/".format(self.gatk_folder)
        
        outfile=""

        if outdir:
            outfile= "{0}/".format(outdir)

        outfile+= "{0}.vcf".format(outprefix)

        command += "GenomeAnalysisTK.jar -T CombineVariants {0} -R {1} ".format(variants_str, reference)

        if ginterval is not None:
            command += "-L {0}".format(ginterval)
        
        if genotypemergeoption is not None:
            command += "--genotypemergeoption {0} ".format(genotypemergeoption)

        if options:
            for opt in options:
                command += "{0} ".format(opt)
        
        if compress is True:
            command += " | bgzip -c > {0}.gz".format(outfile)
        else :
            command += " -o {0}".format(outfile)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            print("Something went wrong while running CombineVariants.\n"
                  "Command used was: %s" % command)
            raise Exception(exc.output)

        if compress is True:
            return outfile+".gz"
        else:
            return outfile
