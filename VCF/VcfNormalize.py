'''
Created on 31 Mar 2017

@author: ernesto
'''

import pdb
import os
import subprocess
import re
from Utils.RunProgram import RunProgram
from collections import namedtuple

class VcfNormalize(object):
    '''
    Normalize variants in a VCF file
    '''


    def __init__(self, vcf, vt_folder=None, vcflib_folder=None, bgzip_folder=None,
                 gatk_folder=None, bcftools_folder=None):
        '''
        Constructor

        Class variables
        ---------------
        vcf : str, Required
             Path to gzipped vcf file
        vt_folder : str, Optional
                      Path to folder containing the vt binary
        vcflib_folder : str, Optional
                        Path to folder containing the different vcflib binaries
        bgzip_folder : str, Optional
                       Path to folder containing the bgzip binary
        gatk_folder : str, Optional
                      Path to folder containing the GATK jar file
        bcftools_folder : str, Optional
                        Path to folder containing the bcftools binary
        '''

        if os.path.isfile(vcf) is False:
            raise Exception("File does not exist")

        self.vcf = vcf
        self.vt_folder = vt_folder
        self.vcflib_folder = vcflib_folder
        self.bgzip_folder = bgzip_folder
        self.gatk_folder = gatk_folder
        self.bcftools_folder = bcftools_folder

    def run_vtnormalize(self,outprefix,reference,compress=False, verbose=False, outdir=None, n=False):
        '''
        Run vt normalize on a vcf file

        Parameters
        ----------
        outprefix : str, required
              prefix for outputfile
        reference : str, required
              path to Fasta file with reference
        compress : boolean, optional
              bgzip compress the normalized VCF
        outdir : str, optional
            If provided, then put output files in this folder
        n : bool, optional
            warns but does not exit when REF is inconsistent
            with reference sequence for non SNPs. Default=False
        verbose : bool, optional
                  if true, then increase verbosity

        Returns
        -------
        A string with path to normalized file
        '''

        if self.vt_folder is None:
            raise Exception("Provide a vt_folder containing the vt binary")

        Arg = namedtuple('Argument', 'option value')

        if outdir:
            outprefix = "{0}/{1}".format(outdir, outprefix)

        outprefix = outprefix+".norm.vcf"

        args=[Arg('-r',reference)]

        parameters=[self.vcf]
        if n is True:
            parameters.append('-n')

        runner=None
        pipelist=None
        if compress is True:
            outprefix += ".gz"
            compressRunner=RunProgram(path=self.bgzip_folder,program='bgzip',parameters=[ '-c', '>', outprefix])
            pipelist=[compressRunner]
        elif compress is None or compress is False:
            args.append(Arg('-o',outprefix))

        runner=RunProgram(path=self.vt_folder, program='vt normalize', args=args, parameters=parameters, downpipe=pipelist)


        if verbose is True:
             print("Command line for running vt normalize is: {0}".format(runner.cmd_line))

        runner.run_checkoutput()
        
        return outprefix

    def run_bcftoolsnorm(self, outprefix, reference, multiallelics=None, type=None, outdir=None,verbose=False):
        '''
        Run bcftools norm on a vcf file

        Parameters
        ----------
        outprefix : str, required
              prefix for outputfile
        reference : str, required
              path to Fasta file with reference
        multiallelic : str, optional
              Operate on multiallelic variants and either split or merge them.
              Possible values are: 'split'/'merge'
        type: str, optional
              If 'multiallelic' is defined then operate on this type of variant.
              Possible values are: snps|indels|both|any
        outdir : str, optional
            If provided, then put output files in this folder
        verbose : bool, optional
                  if true, then increase verbosity

        Returns
        -------
        A string with path to normalized file
        '''

        if outdir:
            outprefix = "{0}/{1}".format(outdir, outprefix)

        outprefix = outprefix+".norm.vcf.gz"

        Arg = namedtuple('Argument', 'option value')
        
        args=[Arg('-f',reference), Arg('-o',outprefix)]

        if multiallelics == "split":
            if type is None: raise Exception("'multiallelics' option is defined, so please provide a 'type' value")
            args.append(Arg('-m',"\'-{0}\'".format(type)))
        elif multiallelics == "merge":
            if type is None: raise Exception("'multiallelics' option is defined, so please provide a 'type' value")
            args.append(Arg('-m',"\'+{0}\'".format(type)))
        else:
            if multiallelics is not None: raise Exception("'multiallelics' value is not recognized: {0}".format(multiallelics))
            
        parameters=[self.vcf,'-Oz']

        runner=RunProgram(path=self.bcftools_folder, program='bcftools norm', args=args, parameters=parameters)

        if verbose is True:
             print("Command line for running bcftools norm is: {0}".format(runner.cmd_line))

        runner.run_checkoutput()

        return outprefix

    def run_vcfallelicprimitives(self, outprefix, compress=True, outdir=None,
                                 keepinfo=True, keepgeno=True, downstream_pipe=None, verbose=None):
        '''
        Run vcfallelicprimitives on a vcf file

        This program is used to decompose complex variants into a canonical SNP and
        indel representation,generating phased genotypes for available samples.

        Parameters
        ----------

        outprefix : str, required
              prefix for outputfiles
        compress : boolean, optional
              bgzip compress the normalized VCF
        outdir : str, optional
            If provided, then put output files in this folder
        keepinfo : bool, optional. Default=True
            Maintain site and allele-level annotations when decomposing.
            Note that in many cases, such as multisample VCFs, these won't
            be valid post-decomposition.  For biallelic loci in single-sample
            VCFs, they should be usable with caution
        keepgeno : bool, optional. Default=True
            Maintain genotype-level annotations when decomposing.  Similar
            caution should be used for this as for keep-info.
        downstream_pipe : str, optional
            If defined, then pipe the output VCF to other tools. 
            i.e. "~/bin/vt/vt sort - | ~/bin/vt/vt uniq -"
        verbose : bool, optional
            if true, then increase verbosity

        Returns
        -------
        A string with path to decomposed file
        '''

        if outdir: 
            outprefix = "{0}/{1}".format(outdir, outprefix)

        outprefix = outprefix+".aprimitives.vcf"

        params=[self.vcf]

        if keepinfo is True:
            params.append('--keep-info')

        if keepgeno is True:
            params.append('--keep-geno')

        if downstream_pipe is not None:
            params.append("| {0}".format(downstream_pipe))

        runner=None
        pipelist=None
        if compress is True:
            outprefix += ".gz"
            compressRunner=RunProgram(path=self.bgzip_folder,program='bgzip',parameters=[ '-c', '>', outprefix])
            pipelist=[compressRunner]
        elif compress is None or compress is False:
            params.extend(['>',outprefix])

        runner=RunProgram(path=self.vcflib_folder, program='vcfallelicprimitives', parameters=params, downpipe=pipelist)

        if verbose is True:
             print("Command line for running vcfallelicprimitives is: {0}".format(runner.cmd_line))

        runner.run_checkoutput()

        return outprefix

    def run_gatk_VariantsToAllelicPrimitives(self, outprefix, reference,
                                             outdir=None, compress=None, verbose=None):
        '''
        Run GATK VariantsToAllelicPrimitives in order to decompose MNPs
         into more basic/primitive alleles

        Parameters
        ----------

        outprefix : str, required
                    prefix for outputfiles
        reference : str, Required
                     Path to fasta file containing the reference
        outdir : str, optional
                   If provided, then put output files in this folder
        compress : boolean, optional
                   bgzip compress the normalized VCF
        verbose : bool, optional
                  if true, then increase verbosity

        Returns
        -------
        A string with path to decomposed file

        '''

        if self.gatk_folder is None:
            raise Exception("Error. I need that the folder containing the GATK "
                            "jar file is defined!")

        if outdir: 
            outprefix = "{0}/{1}".format(outdir, outprefix)

        outprefix = outprefix+".aprimitives.vcf"

        Arg = namedtuple('Argument', 'option value')

        args=[Arg('-T','VariantsToAllelicPrimitives'), Arg('-R',reference),
              Arg('-V',self.vcf), Arg('-o',outprefix) ]
        
        runner=RunProgram(program='java -jar {0}/GenomeAnalysisTK.jar'.format(self.gatk_folder), args=args)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout,stderr=runner.run_popen()

        lines=stderr.split("\n")
        p = re.compile('#* ERROR')
        for i in lines:
            m = p.match(i)
            if m:
                print("Something went wrong while running GATK VariantsToAllelicPrimitives. This was the error message: {0}".format(stderr))
                raise Exception()

        if compress is True:
            compressRunner=RunProgram(path=self.bgzip_folder,program='bgzip',parameters=[ '-c', outprefix, '>', outprefix+".gz"])
            compressRunner.run_checkoutput()
            #delete tmp files
            os.remove(outprefix)
            os.remove(outprefix+".idx")
            outprefix += ".gz"
        
        return outprefix    
