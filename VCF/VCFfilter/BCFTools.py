'''
Created on 27 Feb 2017

@author: ernesto
'''

import os
import glob
import ast
import pdb
from Utils.RunProgram import RunProgram
from collections import namedtuple

class BCFTools(object):
    """
    Class to perform different filtering/selecting operations using BCFTools
    https://samtools.github.io/bcftools/bcftools.html

    Class variables
    ---------------
    bcftools_folder : str, Optional
                      Path to folder containing the bcftools binary
    tabix_folder : str, Optional
                   Path to folder containing the tabix binary
     arg : namedtuple
          Containing a particular argument and its value
    """
    bcftools_folder, tabix_folder= None,None
    arg = namedtuple('Argument', 'option value')

    def __init__(self, vcf: str)->None:
        """
        Constructor

        Parameters
        ----------
        vcf : str
              Path to gzipped vcf file.
        """

        if os.path.isfile(vcf) is False:
            raise Exception("File does not exist")
        self.vcf = vcf

    def subset_vcf(self, prefix: str, bed: str=None, outdir=None,
                   action='exclude', threads: int=1, verbose: bool=False, **kwargs):
        """
        Subset the vcf file using a BED file/region having the coordinates of the
        variants to exclude/include

        Parameters
        ----------
        prefix : str
                 Prefix for outputfiles.
        bed : str, optional
              BED file with coordinates to exclude/include.
        outdir : str, optional
                 If provided, then put output files in this folder.
        action : str, default=exclude
                Exclude or include variants from the bed file passed through the
                bed option.
        threads : int
                  Number of extra output compression threads. Default=1.
        verbose : bool
                  Increase verbosity.
        **kwargs: Arbitrary keyword arguments. Check the `GATK` help for
                  more information.

        Returns
        -------
        prefix : str
                 Path to gzipped VCF file that will have the desired variants excluded/included.
        """
        pdb.set_trace()
        if action != 'include' and action != 'exclude':
            raise Exception("action argument should be either include or exclude")

        if 'r' in kwargs:
            bits = prefix.split(".")
            vcf_ix = bits.index("vcf")

            new = ""
            if 'f' in kwargs:
                new = bits[vcf_ix - 1] + "_" + kwargs['r']+".filt"
            else:
                new = bits[vcf_ix - 1] + "_" + kwargs['r']
            bits[vcf_ix - 1] = new
            prefix = ".".join(bits)

        if outdir:
            prefix = "%s/%s" % (outdir, prefix)

        args = []
        if bed:
            if action == 'exclude':
                args.append(BCFTools.arg('-T', '^{0}'.format(bed)))
            elif action == 'include':
                args.append(BCFTools.arg('-T', '{0}'.format(bed)))
        elif kwargs['r']:
            if action == 'exclude':
                args.append(BCFTools.arg('-t', '^{0}'.format(region)))
            elif action == 'include':
                args.append(BCFTools.arg('-r', '{0}'.format(region)))

         # add now the **kwargs 
        allowed_keys = ['O', 'f', 'r'] # allowed arbitrary args
        args.extend([GATK.arg(f"-{k}", v) for k, v in kwargs.items() if k in allowed_keys])

        args.extend([BCFTools.arg('-o', outprefix), Arg('--threads', threads)])

        cmd = f"{BCFTools.bcftools_folder}/bcftools view" if BCFTools.bcftools_folder else "bcftools view"
        runner = RunProgram(program=cmd,
                            args=args, parameters=[self.vcf])

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout = runner.run_checkoutput()

        return prefix

    def filter(self, name, expression, verbose=None):
        """
        Run bcftools filter on a VCF file

        Parameters
        ----------
        name : str
                 annotate FILTER column with <str>.
        expression : str
                   exclude sites for which expression is true. i.e. 'INFO/DP>24304 | MQ<34'.
        verbose : bool, optional
                  Increase verbosity.

        Returns
        -------
        outfile : str
                Path to the filtered VCF file.
        """

        outfile = self.vcf + ".filtered.vcf.gz"

        Arg = namedtuple('Argument', 'option value')

        args = [Arg('-s', name), Arg('-e', '\'{0}\''.format(expression)), Arg('-o', outfile),
                Arg('-O', 'z')]

        runner = RunProgram(path=self.bcftools_folder, program='bcftools filter',
                            args=args, parameters=[self.vcf])

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout = runner.run_checkoutput()

        return outfile

    def filter_by_variant_type(self, outprefix, v_type="snps", compress=True, biallelic=False,
                               action="select", verbose=None):
        """
        Method to filter a VCF file by variant type. For example, to extract
        only the SNPs

        Parameters
        ----------
        v_type : {'snps','indels','mnps','other','both'}, default=snps
                 Extract/Filter (depending on the value of the 'action'
                 argument) a certain variant type.
        compress : bool, default=True
                   If True then generate a vcf.gz file.
        biallelic : bool, default=False
                    Select only biallelic variants.
        action : {'select', 'exclude'}, default='select'
        outprefix : str
                    Prefix used for the output files.
        verbose : bool, optional
                  Increase verbosity.

        Returns
        -------
        outprefix : str
                 Path to the filtered VCF.
        """

        if v_type not in ('snps', 'indels', 'mnps', 'other', 'both'):
            raise Exception("type value is not valid. Valid values are 'snps'/"
                            "'indels'/'mnps'/'other'/'both'")
        if action not in ('select', 'exclude'):
            raise Exception("action value is not valid. Valid values are 'select' or 'exclude'")

        Arg = namedtuple('Argument', 'option value')

        args = []
        params = []

        if action == "select":
            if v_type != 'both':
                outprefix = outprefix + ".{0}.".format(v_type)
                args.append(Arg('-v', v_type))
        elif action == "exclude":
            if v_type != 'both':
                outprefix = outprefix + ".no{0}.".format(v_type)
                args.append(Arg('-V', v_type))
        if biallelic is True:
            outprefix += "biallelic."
            params.extend(['-m2', '-M2'])

        if compress is True:
            outprefix += "vcf.gz"
            args.extend([Arg('-o', outprefix), Arg('-O', 'z')])
            params.append(self.vcf)
        elif compress is False:
            outprefix += "vcf"
            args.extend([Arg('-o', outprefix), Arg('-O', 'v')])
            params.append(self.vcf)
        elif compress is None:
            raise Exception("'compress' parameter can't be None")

        runner = RunProgram(path=self.bcftools_folder,
                            program='bcftools view',
                            args=args,
                            parameters=params)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout = runner.run_checkoutput()

        return outprefix

    def select_variants(self, outprefix, uncalled=None, threads=1, verbose=None):
        """
        Run bcftools view to select only the variants (exclude the 0|0 genotypes)

        Parameters
        ----------
        outprefix : str
                    Prefix used for the output file.
        uncalled : {'exclude','include'}, optional.
                   Select/Exclude sites with an uncalled genotype.
        threads: int, default=0
                 Number of output compression threads to use in addition to main thread.
        verbose : bool, optional
                  Increase verbosity.

        Returns
        -------
        outfile : str
                Returns the path for the VCF with the selected variants.
        """
        outfile = outprefix + ".onlyvariants.vcf.gz"

        Arg = namedtuple('Argument', 'option value')

        args = [Arg('-o', outfile), Arg('-O', 'z'), Arg('--threads', threads)]

        params = []
        if uncalled == 'exclude':
            params.append('-U')
        elif uncalled == 'include':
            params.append('-u')

        params.append(self.vcf)

        runner = RunProgram(path=self.bcftools_folder,
                            program='bcftools view',
                            args=args,
                            parameters=params)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout = runner.run_checkoutput()

        return outfile
