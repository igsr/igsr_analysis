'''
Created on 27 Feb 2017

@author: ernesto
'''

import os
import subprocess
import json
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
                   create_index=False, verbose=None, action='exclude', apply_filters=None, threads=0):
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

        command = ""
        if self.bcftools_folder:
            command += self.bcftools_folder + "/"

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

        if bed:
            if action == 'exclude':
                command += "bcftools view -T ^{0} ".format(bed)
            elif action == 'include':
                command += "bcftools view -T {0} ".format(bed)
        elif region:
            if action == 'exclude':
                command += "bcftools view -t ^{0} ".format(region)
            elif action == 'include':
                command += "bcftools view -r {0} ".format(region)

        command += "-o {0} -O z {1} --threads {2}".format(outprefix, self.vcf, threads)
        
        if apply_filters is not None:
            command += " -f \"{0}\"".format(apply_filters)

        if verbose is True:
            print("Command is: %s" % command)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            raise Exception(exc.output)

        if create_index is True:
            command_ix = ""
            if self.tabix_folder:
                command_ix += self.tabix_folder + "/"
            command_ix += "tabix {0}".format(outprefix)

            try:
                print("executing cmd: {0}".format(command_ix))
                subprocess.check_output(command_ix, shell=True)
            except subprocess.CalledProcessError as exc:
                raise Exception(exc.output)

        return outprefix

    def filter(self, name, expression):
        '''
        Run bcftools filter on a Vcf file

        Parameters
        ---------
        name :str, Required
                    annotate FILTER column with <str>
        expression :str, Required
                   exclude sites for which expression is true. i.e. 'INFO/DP>24304 | MQ<34'

        Returns
        -------
        Returns the path for the filtered file
        '''
        outfile = self.vcf + ".filtered.vcf.gz"

        command = "{0}/bcftools filter -s{1} -e'{2}' {3} -o {4} -O z".format(
            self.bcftools_folder, name, expression, self.vcf, outfile)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            print("Something went wrong while running Bcftools filter\n"
                  "Command used was: %s" % command)
            raise Exception(exc.output)

        return outfile

    def filter_by_variant_type(self, outprefix, v_type="snps", compress=True, biallelic=False, action="select"):
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
        '''

        if v_type != "snps" and v_type != "indels" and v_type != "mnps" and v_type != "other":
            raise Exception("type value is not valid. Valid values are 'snps'/"
                            "'indels'/'mnps'/'other'")
        if action != "select" and action != "exclude":
            raise Exception("action value is not valid. Valid values are 'select' or 'exclude'")

        command = "{0}/bcftools view ".format(self.bcftools_folder)

        if action == "select":
            outprefix = outprefix + ".{0}.".format(v_type)
            command += "-v {0} ".format(v_type)
        elif action == "exclude":
            outprefix = outprefix + ".no{0}.".format(v_type)
            command += "-V {0} ".format(v_type)
        if biallelic is True:
            outprefix += "biallelic."
            command += "-m2 -M2 "
        
        if compress is True:
            outprefix += "vcf.gz"
            command += "{0} -o {1} -O z".format(self.vcf, outprefix)
        elif compress is False:
            outprefix += "vcf"
            command += "{0} -o {1} -O v".format(self.vcf, outprefix)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            print("Something went wrong while trying to select variants\n")
            print("Command used was: {0}".format(command))
            raise Exception(exc.output)

        return outprefix

    def select_variants(self, outprefix):
        '''
        Run bcftools view to select only the variants (exclude the 0|0 genotypes)

        Parameters
        ---------
        outprefix : str, Required. Prefix used for the output file

        Returns
        -------
        Returns the path for the VCF with only the variants
        '''
        outfile = outprefix + ".onlyvariants.vcf.gz"

        command = "{0}/bcftools view -c1 {1} -o {2} -O z".format(
            self.bcftools_folder, self.vcf, outfile)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            print("Something went wrong while running Bcftools view\n"
                  "Command used was: %s" % command)
            raise Exception(exc.output)

        return outfile
