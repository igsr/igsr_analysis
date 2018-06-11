'''
Created on 27 Feb 2017

@author: ernesto
'''

import os
import subprocess
import json
from Utils.RunProgram import RunProgram
from collections import namedtuple
import glob
import ast
import pdb

class BCFTools(object):
    '''
    Class to perform different filtering/selecting operations using BCFTools
    https://samtools.github.io/bcftools/bcftools.html
    '''
    def __init__(self, vcf, bcftools_folder=None, tabix_folder=None):
        '''
        Constructor

        Class variables
        ---------------
        vcf : str, Required
             Path to gzipped vcf file
        bcftools_folder : str, Required
                         Path to folder containing the bcftools binary
        tabix_folder : str, Optional
                        Path to folder containing the tabix binary
        '''

        if os.path.isfile(vcf) is False:
            raise Exception("File does not exist")

        self.vcf = vcf
        self.bcftools_folder = bcftools_folder
        self.tabix_folder = tabix_folder

    def subset_vcf(self, outprefix, bed=None, region=None, outdir=None,
                   create_index=False, verbose=None, action='exclude', apply_filters=None, threads=1):
        '''
        Subset the vcf file using a BED file/region having the coordinates of the
        variants to exclude/include

        Parameters
        ----------
        bed : string, Optional
              BED file with coordinates to exclude/include
        region : string, Optional
              String with region to consider: chr1, chr1:1000-1500, etc...
        outprefix : str, required
              prefix for outputfiles
        outdir : str, optional
            If provided, then put output files in this folder
        create_index : Boolean, optional
            Generate a tabix index. Default=False
        verbose : Boolean, optional
            verbose
        action: str, optional
            Exclude or include variants from the bed file passed through the
            bed option. Default= exclude
        apply_filters: str, optional
            Apply a filter string: i.e. "PASS,."
        threads: int, optional
            Number of output compression threads to use in addition to main thread. Default=0

        Returns
        -------
        Path to gzipped VCF file that will have the desired variants excluded/included
        '''

        if action != 'include' and action != 'exclude':
            raise Exception("action argument should be either include or exclude")

        if region:
            bits = outprefix.split(".")
            vcf_ix = bits.index("vcf")

            new=""
            if apply_filters is not None:
                new = bits[vcf_ix - 1] + "_" + region+".filt"
            else:
                new = bits[vcf_ix - 1] + "_" + region
            bits[vcf_ix - 1] = new
            outprefix = ".".join(bits)

        if outdir:
            outprefix = "%s/%s" % (outdir, outprefix)

        Arg = namedtuple('Argument', 'option value')

        args=[]
        if bed:
            if action == 'exclude':
                args.append(Arg('-T','^{0}'.format(bed)))
            elif action == 'include':
                args.append(Arg('-T','{0}'.format(bed)))
        elif region:
            if action == 'exclude':
                args.append(Arg('-t','^{0}'.format(region)))
            elif action == 'include':
                args.append(Arg('-r','{0}'.format(region)))

        args.extend([Arg('-o',outprefix), Arg('-O','z'), Arg('--threads',threads)])

        if apply_filters is not None:
            args.append(Arg('-f','\"{0}\"'.format(apply_filters)))

        runner=RunProgram(path=self.bcftools_folder, program='bcftools view', args=args, parameters=[self.vcf])

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout=runner.run_checkoutput()

        return outprefix

    def filter(self, name, expression, verbose=None):
        '''
        Run bcftools filter on a Vcf file

        Parameters
        ---------
        name :str, Required
                    annotate FILTER column with <str>
        expression :str, Required
                   exclude sites for which expression is true. i.e. 'INFO/DP>24304 | MQ<34'
        verbose : Boolean, optional
                  Increase verbosity

        Returns
        -------
        Returns the path for the filtered file
        '''

        outfile = self.vcf + ".filtered.vcf.gz"

        Arg = namedtuple('Argument', 'option value')

        args=[Arg('-s',name),Arg('-e',expression),Arg('-o',outfile),Arg('-O','z')]

        runner=RunProgram(path=self.bcftools_folder, program='bcftools filter', args=args, parameters=[self.vcf])

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout=runner.run_checkoutput()

        return outfile

    def filter_by_variant_type(self, outprefix, v_type="snps", compress=True, biallelic=False, action="select", verbose=None):
        '''
        Method to filter a VCF file by variant type. For example, to extract only the SNPs
        
        Parameters
       ----------
        v_type : str, Required. Valid values are 'snps'/'indels','mnps','other'. Default=snps
                 Extract/Filter (depending on the value of the 'action'
                 argument) a certain variant type
        compress : bool, Optional
                   If True then generate a vcf.gz file. Default=True 
        biallelic : bool, Optional
                    Select only biallelic variants. Default=False
        action : str, Required. Valid values are 'select'/'exclude'. Default=select
        outprefix : str, Required. Prefix used for the output files
        verbose : Boolean, optional
                  Increase verbosity
        '''

        if v_type != "snps" and v_type != "indels" and v_type != "mnps" and v_type != "other":
            raise Exception("type value is not valid. Valid values are 'snps'/"
                            "'indels'/'mnps'/'other'")
        if action != "select" and action != "exclude":
            raise Exception("action value is not valid. Valid values are 'select' or 'exclude'")

        Arg = namedtuple('Argument', 'option value')

        args=[]
        params=[]

        if action == "select":
            outprefix = outprefix + ".{0}.".format(v_type)
            args.append(Arg('-v',v_type))
        elif action == "exclude":
            outprefix = outprefix + ".no{0}.".format(v_type)
            args.append(Arg('-V',v_type))
        if biallelic is True:
            outprefix += "biallelic."
            params.extend(['-m2','-M2'])
        
        if compress is True:
            outprefix += "vcf.gz"
            args.extend([Arg('-o', outprefix),Arg('-O','z')])
            params.append(self.vcf)
        elif compress is False:
            outprefix += "vcf"
            args.extend([Arg('-o', outprefix),Arg('-O','v')])
            params.append(self.vcf)

        runner=RunProgram(path=self.bcftools_folder, program='bcftools view', args=args, parameters=params)
        
        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout=runner.run_checkoutput()

        return outprefix

    def select_variants(self, outprefix, uncalled=None, threads=1, verbose=None):
        '''
        Run bcftools view to select only the variants (exclude the 0|0 genotypes)

        Parameters
        ---------
        outprefix : str, Required. Prefix used for the output file
        uncalled : str, optional. Select/Exclude sites with an uncalled genotype. 
                   Possible values are: 'exclude', 'include'
        threads: int, optional
                 Number of output compression threads to use in addition to main thread. Default=0
        verbose : Boolean, optional
                  Increase verbosity

        Returns
        -------
        Returns the path for the VCF with only the variants
        '''
        outfile = outprefix + ".onlyvariants.vcf.gz"

        Arg = namedtuple('Argument', 'option value')

        args=[Arg('-o',outfile),Arg('-O','z'),Arg('--threads',threads)]

        params=[]
        if uncalled=='exclude':
            params.append('-U')
        elif uncalled=='include':
            params.append('-u')

        params.append(self.vcf)

        runner=RunProgram(path=self.bcftools_folder, program='bcftools view', args=args, parameters=params)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout=runner.run_checkoutput()

        return outfile
