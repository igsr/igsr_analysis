"""VcfQC module.

This module is used to do a series of VCF QC tasks on a VCF file.

@author: Ernesto Lowy
"""

import os
import subprocess
import tempfile
import pdb

from Utils.RunProgram import RunProgram
from collections import namedtuple

class VcfQC(object):
    '''
    Class to do quality assessment on a VCF file
    '''


    def __init__(self, vcf, bgzip_folder=None, bcftools_folder=None,
                 bedtools_folder=None, picard_folder=None, r_folder=None,
                 r_scripts=None, tabix_folder=None):
        '''
        Constructor

        Parameters
        ----------
        vcf : str
              Path to gzipped vcf file.
        bgzip_folder : str, optional
                       Path to folder containing the bgzip binary.
        bcftools_folder : str, optional
                          Path to folder containing the bcftools binary.
        bedtools_folder : str, optional
                          Path to the folder containing the bedtools binary.
        picard_folder : str, optional
                        Path to folder containing the picard binary.
        r_folder : str, optional
                   Path to folder containing the R binary.
        r_scripts : str, optional
                    Path to folder containing the R scripts required for
                    constructing some of the plots used by this class (i.e. plot_variants.R).
        tabix_folder : str, optional
                        Path to folder containing the tabix binary.
        '''

        if os.path.isfile(vcf) is False:
            raise IOError("File does not exist")

        self.vcf = vcf
        self.bgzip_folder = bgzip_folder
        self.bcftools_folder = bcftools_folder
        self.bedtools_folder = bedtools_folder
        self.picard_folder = picard_folder
        self.r_folder = r_folder
        self.r_scripts = r_scripts
        self.tabix_folder = tabix_folder

    def calc_concordance(self, truth_vcf, truth_sample, call_sample, outprefix,
                         outdir=None, intervals=None, verbose=None):
        '''
        Method to calculate the genotype concordance between VcfQC.vcf and Truth VCF.
        It will run Picard's GenotypeConcordance

        Parameters
        ----------
        truth_vcf : str
                    The VCF containing the truth sample.
        truth_sample : str
                       The name of the truth sample within the truth VCF.
        call_sample : str
                      The name of the call sample within the call VCF.
        outprefix : str
                    String used as the prefix in the output file.
        outdir : str, optional
                 If provided, then put output files in this folder.
        intervals : str
                    One or more interval list files that will be used to limit the
                    genotype concordance.
        verbose : bool, optional
                  Ff true, then print the command line used for running this program.

        Returns
        -------
        GTPconcordance object
        '''

        if self.picard_folder is None:
            raise Exception("Folder containing Picard jar file is required")

        if outdir:
            outprefix = "%s/%s" % (outdir, outprefix)

        Arg = namedtuple('Argument', 'option value')

        args=[Arg('TRUTH_VCF',truth_vcf), Arg('CALL_VCF',self.vcf), Arg('TRUTH_SAMPLE',truth_sample), 
              Arg('CALL_SAMPLE',call_sample), Arg('O',outprefix)] 

        if intervals:
            args.append(Arg('INTERVALS',intervals))

        runner=RunProgram(program='java -jar {0}/picard.jar GenotypeConcordance'.format(self.picard_folder), args=args, arg_sep='=')

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout=runner.run_checkoutput()
        
        gtp_con = GTPconcordance(summary_metrics_file=outprefix+\
                                 ".genotype_concordance_summary_metrics")

        return gtp_con

    def get_chros(self, filter_str=None, chr_f=None, verbose=None):
        '''
        Method to get a list of chromosomes present in a file

        Parameters
        ----------
        filter_str : str, optional
                     If defined, apply this filter string so bcftools view
                     apply it before fetching the chros.
        chr_f : str, optional
                Path to file with a list of chromosomes (one per line).
                If provided, the chros in the file will be compared with the
                chromosomes in self.vcf.
        verbose : bool, optional
                  Increase verbosity.

        Returns
        -------
        dict
            Dict with a key named 'in_vcf' and whose values are the chros that are present in self.vcf.

            If list_of_chros is defined, then it will also add 3 keys to the dict:
                 'both' whose values will be the chros present in self.vcf and in 'chr_f'
                 'in_A' whose values will be the chros PRESENT in self.vcf and NOT in 'chr_f'
                 'in_B' whose values will be the chros NOT present in self.vcf and PRESENT in 'chr_f'.
        '''

        params=['--no-header',self.vcf, "|cut -f1 |uniq"]
        
        Arg = namedtuple('Argument', 'option value')

        args=[]

        if filter_str != None:
            args.append(Arg('-f', filter_str))

        runner=RunProgram(path=self.bcftools_folder, program='bcftools view', args=args, parameters=params)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        out_str = ""

        out=runner.run_checkoutput()
        out_str = out.decode("utf-8")
        out_str = out_str.rstrip('\n')

        list_of_chros = out_str.split("\n")

        chr_list_f = []
        if chr_f != None:
            #parse file with chros
            chr_file = open(chr_f, 'r')
            chr_list_f = chr_file.read().splitlines()

        both = set(list_of_chros) & set(chr_list_f)
        in_a = set(list_of_chros) - set(chr_list_f)
        in_b = set(chr_list_f) - set(list_of_chros)

        return {
            'in_vcf' : list_of_chros,
            'both' : list(both),
            'in_A' : list(in_a),
            'in_B' : list(in_b)
            }

    def number_variants_in_region(self, region, outprefix, verbose=None):
        '''
        Method to get the number of variants in a particular region/s

        Parameters
        ----------
        region : str
                 String with path to the BED file containing the regions for
                 which the number will be calculated.
        outprefix : str
                    Prefix for outfile.
        verbose : bool, optional
                  Increase verbosity.

        Returns
        -------
        filename
                File with the number of variants for each particular region.
        '''

        outprefix = outprefix+".counts"
        
        params=['-counts','>',outprefix]

        Arg = namedtuple('Argument', 'option value')

        args=[Arg('-a',region), Arg('-b', self.vcf)]

        runner=RunProgram(path=self.bedtools_folder, program='bedtools coverage', args=args, parameters=params)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout=runner.run_checkoutput()

        return outprefix

    def run_CollectVariantCallingMetrics(self, outprefix, truth_vcf, intervals=None, verbose=None):

        '''
        Method to run Picard's CollectVariantCallingMetrics on a VcfQC object.

        Parameters
        ----------
        outprefix : str
                    Prefix for outfiles: prefix.variant_calling_detail_metrics
                    and prefix.variant_calling_summary_metrics.
        truth_vcf : str
                    Reference VCF file.
        intervals : str, optional
                    Target intervals to restrict analysis to.
        verbose : bool, optional
                  Increase verbosity.

        Returns
        -------
        CollectVCallingMetrics object
        '''

        if self.picard_folder is None:
            raise Exception("Provide a picard folder") 


        Arg = namedtuple('Argument', 'option value')

        args=[Arg('I',self.vcf), Arg('O',outprefix), Arg('DBSNP',truth_vcf)]

        if intervals:
            args.append(Arg('TI', intervals))

        runner=RunProgram(program='java -jar {0}/picard.jar CollectVariantCallingMetrics'.format(self.picard_folder), args=args, arg_sep="=")
            
        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout=runner.run_checkoutput()

        #create CollectVCallingMetrics object with the output files
        cvcmetrics = CollectVCallingMetrics(vc_detail_metrics_file=outprefix\
                                            +".variant_calling_detail_metrics",
                                            vc_summary_metrics_file=outprefix\
                                            +".variant_calling_summary_metrics")

        return cvcmetrics
    
    def plot_variant_density(self,length,genome,outprefix,plot_params):
        '''
        Method used to plot the variant density along the genome for this VCF file

        This method will calculate the variant density (SNPs, INDELs or both)
        for windows of a certain length (in bp) defined by the user and will
        plot this density along the chromosomes

        Parameters
        ----------
        length : int
                 Length in bp of the genomic windows for which the number
                 of variants will be calculated.
        genome : str
                 File with genome sizes.
        outprefix : str
                    Prefix for output file (png file with variant density plot).
        plot_params : dict
                      Dictionary containing the graphical parameters used for plot_variants.R:
                      {
                       'height': int,
                       'chromName_cex': double
                       }.

        Returns
        -------
        filename
                Returns png file with the density plot
        '''

        #run bedtools makewindows to slice the genome
        windows_tmp = tempfile.NamedTemporaryFile()

        command1 = ""

        if self.bedtools_folder:
            command1 += self.bedtools_folder+"/"

        command1 += "bedtools makewindows -g {0} -w {1} > {2}"\
        .format(genome, length, windows_tmp.name)

        try:
            subprocess.check_output(command1, shell=True)
        except subprocess.CalledProcessError as exc:
            print("Command used was: {0}".format(command1))
            raise Exception(exc.output)

        #run bedtools coverage
        cov_tmp = tempfile.NamedTemporaryFile()

        command2 = ""

        if self.bedtools_folder:
            command2 += self.bedtools_folder+"/"

        command2 += "bedtools coverage -a {0} -b {1} -counts > {2}"\
        .format(windows_tmp.name, self.vcf, cov_tmp.name)

        try:
            subprocess.check_output(command2, shell=True)
        except subprocess.CalledProcessError as exc:
            raise Exception(exc.output)

        #now, use an R script to calculate the density and plot it along the chromosomes
        command3 = ""

        if self.r_folder:
            command3 += self.r_folder+"/"

        if self.r_scripts is None:
            raise Exception("Error! I need that the r_scripts arg is defined.")
        
        command3 +="R --slave < {0}/plot_variants.R --args {1} {2} {3} {4} {5}".format(self.r_scripts, cov_tmp.name,
                                                                                       length, outprefix,plot_params['height'],
                                                                                       plot_params['chromName_cex'])
        try:
            subprocess.check_output(command3, shell=True)
        except subprocess.CalledProcessError as exc:
            raise Exception(exc.output)

        # Clean up the files
        windows_tmp.close()
        cov_tmp.close()

        return outprefix+".png"

    def stats(self, outpath, filter_str=None, region=None, region_file=None, verbose=None):
        '''
        Run bcftools stats on the VCF file

        Parameters
        ----------
        outpath : str
                  output path
        filter_str : str, optional. 
                     Example:  PASS,.
                     Apply filters when calculating the stats.
        region : str, optional
                 Example: chr20
                 Region used for calculating the stats.
        region_file : filename, optional
                      BED file with the regions that will be analyzed.
        verbose : bool, optional

        Returns
        -------
        BcftoolsStats object
        '''

        Arg = namedtuple('Argument', 'option value')
          
        args=[]

        if region != None:
            outpath = "{0}.{1}".format(outpath, region)
            args.append(Arg('-r',region))

        if region_file != None:
            args.append(Arg('-R',region_file))

        if filter_str != None:
            outpath = outpath+".filter_str"
            args.append(Arg('-f',filter_str))

        outpath = outpath+".stats"
        
        params=[self.vcf,'>',outpath]

        runner=RunProgram(path=self.bcftools_folder, program='bcftools stats', args=args, parameters=params)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        runner.run_checkoutput()

        stats = BcftoolsStats(filename=outpath)

        with open(outpath) as fi:
            d = {}
            for line in fi:
                line = line.rstrip('\n')
                if line.startswith('SN\t'):
                    key = line.split('\t')[2]
                    value = int(line.split('\t')[3])
                    d[key] = value
                elif line.startswith('TSTV\t'):
                    ts_tv = line.split('\t')[4]
                    ts_tv_1stalt = line.split('\t')[7]
                    stats.ts_tv = ts_tv
                    stats.ts_tv_1stalt = ts_tv_1stalt
                elif line.startswith('SiS\t'):
                    no_singleton_snps = line.split('\t')[3]
                    stats.no_singleton_snps = no_singleton_snps
                    
            stats.summary_numbers = d
        return stats


class BcftoolsStats(object):
    '''
    Class to store the results of running BCFtools stats on a VCF file
    '''

    def __init__(self, filename=None, summary_numbers=None, ts_tv=None,
                 ts_tv_1stalt=None, no_singleton_snps=None):
        '''
        Constructor

        Parameters
        ----------
        filename : str
                   Filename of the VCF that was used to run bcftools stats
        summary_numbers : dict
                          Dictionary containing the basic stats. i.e.:
                              number of samples:      1
                              number of records:      1867316
                              .....
        ts_tv : float
                ts/tv ratio
        ts_tv_1stalt : float
                       ts/tv (1st ALT)
        no_singleton_snps : int
        '''

        self.filename = filename
        self.summary_numbers = summary_numbers
        self.ts_tv = ts_tv
        self.ts_tv_1stalt = ts_tv_1stalt
        self.no_singleton_snps = no_singleton_snps

    def __str__(self):
        sb = []
        for key in self.__dict__:
            sb.append("{key}='{value}'".format(key=key, value=self.__dict__[key]))

        return ', '.join(sb)

    def __repr__(self):
        return self.__str__()

class GTPconcordance(object):
    '''
    Class to store the results of running picard's GenotypeConcordance on a VCF file
    '''

    def __init__(self, summary_metrics_file, summary_metrics_snps=None,
                 summary_metrics_indels=None):
        '''
        Constructor

        Parameters
        ----------
        summary_metrics_file : str
                               Filepath to *.genotype_concordance_summary_metrics generated by
                               Picard's GenotypeConcordance
        summary_metrics_snps : dict, optional
                               Dict with the summary metrics found in the
                               *.genotype_concordance_summary_metrics file (only for SNPs)
        summary_metrics_indels : dict, optional
                                 Dict with the summary metrics found in the
                                 *.genotype_concordance_summary_metrics file (only for Indels)
        '''

        sm_snps = {}
        sm_indels = {}
        with open(summary_metrics_file) as f:
            keys = []
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('VARIANT_TYPE\t'):
                    keys = line.split('\t')
                elif line.startswith('SNP\t'):
                    values = line.split('\t')
                    sm_snps = dict(zip(keys, values))
                elif line.startswith('INDEL\t'):
                    values = line.split('\t')
                    sm_indels = dict(zip(keys, values))

        self.summary_metrics_snps = sm_snps
        self.summary_metrics_indels = sm_indels

    def __str__(self):
        sb = []
        for key in self.__dict__:
            sb.append("{key}='{value}'".format(key=key, value=self.__dict__[key]))

        return ', '.join(sb)

    def __repr__(self):
        return self.__str__()

class CollectVCallingMetrics(object):
    '''
    Class to store the results of running picard's CollectVariantCallingMetrics on a VCF file
    '''

    def __parse_CollectVariantCallingMetrics_file(self, file):
        '''
        Private function to parse output of CollectVariantCallingMetrics

        Returns
        -------
        dict
            A dictionary with the information in the file
        '''
        out_dict = {}
        with open(file) as f:
            keys = []
            seen = False
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('TOTAL_SNPS\t') or line.startswith('SAMPLE_ALIAS\t'):
                    seen = True
                    keys = line.split('\t')
                elif seen is True:
                    values = line.split('\t')
                    out_dict = dict(zip(keys, values))
                    break
        if seen is False:
            raise Exception("Error parsing {0}. No lines were found in this file! ".format(file))
        return out_dict

    def __init__(self, vc_detail_metrics_file, vc_summary_metrics_file,
                 vc_detail_metrics=None, vc_summary_metrics=None):
        '''
        Constructor

        Parameters
        ----------
        vc_detail_metrics_file : str
                                 Filepath to *.variant_calling_detail_metrics
                                 generated by Picard's CollectVariantCallingMetrics
        vc_summary_metrics_file : str
                                  Filepath to *.variant_calling_summary_metrics
                                  generated by Picard's CollectVariantCallingMetrics
        vc_detail_metrics : dict, optional
                            Dictionary with information in file
                            *.variant_calling_detail_metrics
        vc_summary_metrics : dict, optional
                             Dictionary with information in file
                             *.variant_calling_summary_metrics

        '''

        if os.path.isfile(vc_detail_metrics_file) is False:
            raise Exception("File *.variant_calling_detail_metrics does not exist")

        if os.path.isfile(vc_summary_metrics_file) is False:
            raise Exception("File *.variant_calling_summary_metrics does not exist")

        vc_detail_metrics = self.__parse_CollectVariantCallingMetrics_file(vc_detail_metrics_file)
        vc_summary_metrics = self.__parse_CollectVariantCallingMetrics_file(vc_summary_metrics_file)

        self.vc_detail_metrics_file = vc_detail_metrics_file
        self.vc_summary_metrics_file = vc_summary_metrics_file
        self.vc_detail_metrics = vc_detail_metrics
        self.vc_summary_metrics = vc_summary_metrics

    def __str__(self):
        sb = []
        for key in self.__dict__:
            sb.append("{key}='{value}'".format(key=key, value=self.__dict__[key]))

        return ', '.join(sb)

    def __repr__(self):
        return self.__str__()
