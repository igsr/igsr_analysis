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
        vcf : Path to gzipped vcf file.
        """

        if os.path.isfile(vcf) is False:
            raise Exception("File does not exist")
        self.vcf = vcf

    def subset_vcf(self, prefix: str, bed: str=None,
                   action: str='exclude', threads: int=1, verbose: bool=False, **kwargs)->str:
        """
        Subset the vcf file using a BED file/region having the coordinates of the
        variants to exclude/include

        Parameters
        ----------
        prefix : Prefix for outputfiles.
        bed : BED file with coordinates to exclude/include.
        action : Exclude or include variants from the bed file passed through the
                 bed option.
        threads : Number of extra output compression threads. Default=1.
        verbose : Increase verbosity.
        **kwargs: Arbitrary keyword arguments. Check the `bcftools view` help for
                  more information.

        Returns
        -------
        prefix : Path to gzipped VCF file that will have the desired variants excluded/included.
        """
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

        args = []
        if bed:
            if action == 'exclude':
                args.append(BCFTools.arg('-T', '^{0}'.format(bed)))
            elif action == 'include':
                args.append(BCFTools.arg('-T', '{0}'.format(bed)))
        elif kwargs['r']:
            if action == 'exclude':
                args.append(BCFTools.arg('-t', '^{0}'.format(kwargs['r'])))
            elif action == 'include':
                args.append(BCFTools.arg('-r', '{0}'.format(kwargs['r'])))

         # add now the **kwargs 
        allowed_keys = ['O', 'f'] # allowed arbitrary args
        args.extend([BCFTools.arg(f"-{k}", v) for k, v in kwargs.items() if k in allowed_keys])

        args.extend([BCFTools.arg('-o', prefix), BCFTools.arg('--threads', threads)])

        cmd = f"{BCFTools.bcftools_folder}/bcftools view" if BCFTools.bcftools_folder else "bcftools view"
        runner = RunProgram(program=cmd,
                            args=args, parameters=[self.vcf])

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout = runner.run_checkoutput()

        return prefix

    def filter(self, expression: str, verbose: bool=False, **kwargs)->str:
        """
        Run bcftools filter on a VCF file

        Parameters
        ----------
        verbose : Increase verbosity.
        **kwargs: Arbitrary keyword arguments. Check the `bcftools` help for
                  more information.

        Returns
        -------
        outfile : Path to the filtered VCF file.
        """

        outfile = None
        if 'O' in kwargs:
            if kwargs['O'] == 'z':
                outfile = f"{self.vcf}.filtered.vcf.gz"
        else:
            outfile = f"{self.vcf}.vcf"
        

        arguments = [BCFTools.arg('-o', outfile)]

        # add now the **kwargs        
        allowed_keys = ['s', 'e', 'O'] # allowed arbitrary args
        arguments.extend([BCFTools.arg(f"-{k}", v) for k, v in kwargs.items() if k in allowed_keys])

        cmd = f"{BCFTools.bcftools_folder}/bcftools view" if BCFTools.bcftools_folder else "bcftools view"
        runner = RunProgram(program=cmd,
                            args=arguments, parameters=[self.vcf])

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout = runner.run_checkoutput()

        return outfile

    def filter_by_variant_type(self, prefix: str, v_type="snps", compress: bool=True, biallelic: bool=False,
                               action="select", verbose=False)->str:
        """
        Method to filter a VCF file by variant type. For example, to extract
        only the SNPs

        Parameters
        ----------
        prefix : Prefix for output files.
        v_type : {'snps','indels','mnps','other','both'}, default=snps
                 Extract/Filter (depending on the value of the 'action'
                 argument) a certain variant type.
        compress : If True then generate a vcf.gz file.
        biallelic : Select only biallelic variants.
        action : {'select', 'exclude'}, default='select'
        verbose : Increase verbosity.

        Returns
        -------
        prefix : Path to the filtered VCF.
        """

        if v_type not in ('snps', 'indels', 'mnps', 'other', 'both'):
            raise Exception("type value is not valid. Valid values are 'snps'/"
                            "'indels'/'mnps'/'other'/'both'")
        if action not in ('select', 'exclude'):
            raise Exception("action value is not valid. Valid values are 'select' or 'exclude'")

        args = [] 
        params = []

        if action == "select":
            if v_type != 'both':
                prefix = prefix + ".{0}.".format(v_type)
                args.append(BCFTools.arg('-v', v_type))
        elif action == "exclude":
            if v_type != 'both':
                prefix = outprefix + ".no{0}.".format(v_type)
                args.append(BCFTools.arg('-V', v_type))
        if biallelic is True:
            prefix += "biallelic."
            params.extend(['-m2', '-M2'])

        if compress is True:
            prefix += "vcf.gz"
            args.extend([BCFTools.arg('-o', prefix), BCFTools.arg('-O', 'z')])
            params.append(self.vcf)
        elif compress is False:
            prefix += "vcf"
            args.extend([BCFTools.arg('-o', prefix), BCFTools.arg('-O', 'v')])
            params.append(self.vcf)
        elif compress is None:
            raise Exception("'compress' parameter can't be None")

        cmd = f"{BCFTools.bcftools_folder}/bcftools view" if BCFTools.bcftools_folder else "bcftools view"
        runner = RunProgram(program=cmd,
                            args=args,
                            parameters=params)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout = runner.run_checkoutput()

        return prefix

    def select_variants(self, outprefix : str, uncalled : str = None, threads : int = 1, verbose : bool = False)->str:
        """
        Run bcftools view to select only the variants (exclude the 0|0 genotypes)

        Parameters
        ----------
        outprefix : Prefix used for the output file.
        uncalled : {'exclude','include'},
                   Select/Exclude sites with an uncalled genotype.
        threads: Number of output compression threads to use in addition to main thread.
        verbose : Increase verbosity.

        Returns
        -------
        outfile : Returns the path for the VCF with the selected variants.
        """
        outfile = outprefix + ".onlyvariants.vcf.gz"

        args = [BCFTools.arg('-o', outfile), BCFTools.arg('-O', 'z'), BCFTools.arg('--threads', threads)]

        params = []
        if uncalled == 'exclude':
            params.append('-U')
        elif uncalled == 'include':
            params.append('-u')

        params.append(self.vcf)

        cmd = f"{BCFTools.bcftools_folder}/bcftools view" if BCFTools.bcftools_folder else "bcftools view"
        runner = RunProgram(program=cmd,
                            args=args,
                            parameters=params)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout = runner.run_checkoutput()

        return outfile
