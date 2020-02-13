"""VcfNormalize module.

This module is used to normalize a file in the VCF format
(see https://samtools.github.io/hts-specs/VCFv4.2.pdf).
Normalizing a file is basically used to make different VCFs
comparable, and it is highly recommended when you are
benchmarcking your call set using a reference call set.
Normalization is specially important for INDELs, as the
same INDEL could be represented in different ways depending
on the caller.
"""

import os
from collections import namedtuple
from Utils.RunProgram import RunProgram

class VcfNormalize:
    '''
    Normalize variants in a VCF file
    '''

    def __init__(self, vcf, vt_folder=None, vcflib_folder=None, bgzip_folder=None,
                 gatk_folder=None, bcftools_folder=None):
        '''
        Constructor

        Parameters
        ----------
        vcf : str
              Path to gzipped vcf file
        vt_folder : str, optional
                    Path to folder containing the vt binary.
        vcflib_folder : str, optional
                        Path to folder containing the different vcflib binaries.
        bgzip_folder : str, optional
                       Path to folder containing the bgzip binary.
        gatk_folder : str, optional
                      Path to folder containing the GATK jar file.
        bcftools_folder : str, optional
                          Path to folder containing the bcftools binary.
        '''

        if os.path.isfile(vcf) is False:
            raise Exception("File does not exist")

        self.vcf = vcf
        self.vt_folder = vt_folder
        self.vcflib_folder = vcflib_folder
        self.bgzip_folder = bgzip_folder
        self.gatk_folder = gatk_folder
        self.bcftools_folder = bcftools_folder

    def run_vtnormalize(self, outprefix, reference, compress=False,
                        verbose=False, outdir=None, n=False):
        '''
        Run vt normalize on a vcf file

        Parameters
        ----------
        outprefix : str
                    Prefix for outputfile.
        reference : str
                    Path to Fasta file with reference.
        compress : boolean, optional
                   bgzip compress the normalized VCF.
        outdir : str, optional
                 If provided, then put output files in this folder.
        n : bool, optional
            warns but does not exit when REF is inconsistent
            with reference sequence for non SNPs. Default=False.
        verbose : bool, optional
                  if true, then increase verbosity.

        Returns
        -------
        str
           A string with path to normalized file.
        '''

        if self.vt_folder is None:
            raise Exception("Provide a vt_folder containing the vt binary")

        Arg = namedtuple('Argument', 'option value')

        if outdir:
            outprefix = "{0}/{1}".format(outdir, outprefix)

        outprefix = outprefix+".norm.vcf"

        args = [Arg('-r', reference)]

        parameters = [self.vcf]
        if n is True:
            parameters.append('-n')

        runner = None
        pipelist = None
        if compress is True:
            outprefix += ".gz"
            compressRunner = RunProgram(path=self.bgzip_folder,
                                        program='bgzip',
                                        parameters=['-c', '>',
                                                    outprefix])
            pipelist = [compressRunner]
        elif compress is None or compress is False:
            args.append(Arg('-o', outprefix))

        runner = RunProgram(path=self.vt_folder,
                            program='vt normalize',
                            args=args,
                            parameters=parameters,
                            downpipe=pipelist)

        if verbose is True:
            print("Command line for running vt normalize is: {0}".format(runner.cmd_line))

        runner.run_checkoutput()

        return outprefix

    def run_bcftoolsnorm(self, outprefix, reference, multiallelics=None,
                         type=None, outdir=None, verbose=False):
        '''
        Run bcftools norm on a vcf file

        Parameters
        ----------
        outprefix : str
                    Prefix for outputfile.
        reference : str
                    Path to Fasta file with reference.
        multiallelic : str, optional
                       Operate on multiallelic variants and either split or merge them.
                       Possible values are: 'split'/'merge'
        type: : {'snps', 'indels', 'both', 'any'}, optional
                If 'multiallelic' is defined then operate on this type of variant.
        outdir : str, optional
                 If provided, then put output files in this folder.
        verbose : bool, optional
                  Ff true, then increase verbosity.

        Returns
        -------
        str
           A string with path to normalized file.
        '''

        if outdir:
            outprefix = "{0}/{1}".format(outdir, outprefix)

        outprefix = outprefix+".norm.vcf.gz"

        Arg = namedtuple('Argument', 'option value')

        args = [Arg('-f', reference), Arg('-o', outprefix)]

        if multiallelics == "split":
            if type is None:
                raise Exception("'multiallelics' option is defined, "
                                "so please provide a 'type' value")
            args.append(Arg('-m', "\'-{0}\'".format(type)))
        elif multiallelics == "merge":
            if type is None:
                raise Exception("'multiallelics' option is defined,"
                                " so please provide a 'type' value")
            args.append(Arg('-m', "\'+{0}\'".format(type)))
        else:
            if multiallelics is not None:
                raise Exception("'multiallelics' value is not "
                                "recognized: {0}".format(multiallelics))

        parameters = [self.vcf, '-Oz']

        runner = RunProgram(path=self.bcftools_folder,
                            program='bcftools norm',
                            args=args,
                            parameters=parameters)

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
        outprefix : str
                    Prefix for outputfiles.
        compress : bool, optional
                   Bgzip compress the normalized VCF.
        outdir : str, optional
                 If provided, then put output files in this folder.
        keepinfo : bool, optional
                   Maintain site and allele-level annotations when decomposing.
                   Note that in many cases, such as multisample VCFs, these won't
                   be valid post-decomposition.  For biallelic loci in single-sample
                   VCFs, they should be usable with caution. Default=True.
        keepgeno : bool, optional
                   Maintain genotype-level annotations when decomposing. Similar
                   caution should be used for this as for keep-info. Default=True.
        downstream_pipe : str, optional
                          If defined, then pipe the output VCF to other tools.
                          i.e. "~/bin/vt/vt sort - | ~/bin/vt/vt uniq -".
        verbose : bool, optional
                  if true, then increase verbosity.

        Returns
        -------
        str
           A string with path to decomposed file
        '''

        if outdir:
            outprefix = "{0}/{1}".format(outdir, outprefix)

        outprefix = outprefix+".aprimitives.vcf"

        params = [self.vcf]

        if keepinfo is True:
            params.append('--keep-info')

        if keepgeno is True:
            params.append('--keep-geno')

        if downstream_pipe is not None:
            params.append("| {0}".format(downstream_pipe))

        runner = None
        pipelist = None
        if compress is True:
            outprefix += ".gz"
            compressRunner = RunProgram(path=self.bgzip_folder,
                                        program='bgzip',
                                        parameters=['-c', '>', outprefix])
            pipelist = [compressRunner]
        elif compress is None or compress is False:
            params.extend(['>', outprefix])

        runner = RunProgram(path=self.vcflib_folder,
                            program='vcfallelicprimitives',
                            parameters=params,
                            downpipe=pipelist)

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
        outprefix : str
                    Prefix for outputfiles.
        reference : str
                    Path to fasta file containing the reference.
        outdir : str, optional
                 If provided, then put output files in this folder.
        compress : boolean, optional
                   Bgzip compress the normalized VCF.
        verbose : bool, optional
                  Ff true, then increase verbosity.

        Returns
        -------
        str
           A string with path to decomposed file.
        '''

        if self.gatk_folder is None:
            raise Exception("Error. I need that the folder containing the GATK "
                            "jar file is defined!")

        if outdir:
            outprefix = "{0}/{1}".format(outdir, outprefix)

        outprefix = outprefix+".aprimitives.vcf"

        Arg = namedtuple('Argument', 'option value')

        args = [Arg('-T', 'VariantsToAllelicPrimitives'), Arg('-R', reference),
                Arg('-V', self.vcf), Arg('-o', outprefix)]

        runner = RunProgram(program='java -jar {0}/GenomeAnalysisTK.jar'.
                            format(self.gatk_folder), args=args)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout, stderr, is_error = runner.run_popen()

        if compress is True:
            compressRunner = RunProgram(path=self.bgzip_folder, program='bgzip',
                                        parameters=['-c', outprefix, '>', outprefix+".gz"])
            compressRunner.run_checkoutput()
            #delete tmp files
            os.remove(outprefix)
            os.remove(outprefix+".idx")
            outprefix += ".gz"
        elif compress is False:
            return outprefix
        else:
            raise Exception("'compress' parameter is not valid")

        return outprefix
