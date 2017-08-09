'''
Created on 13 Feb 2017

@author: ernesto
'''

import os
import subprocess
import tempfile

class VcfQC(object):
    '''
    Class to do quality assessment on a VCF file
    '''


    def __init__(self, vcf, bgzip_folder=None, bcftools_folder=None,
                 bedtools_folder=None, picard_folder=None, r_folder=None,
                 r_scripts=None, tabix_folder=None):
        '''
        Constructor

        Class variables
        ---------------
        vcf : str, Required
             Path to gzipped vcf file
        bgzip_folder : str, Optional
                      Path to folder containing the bgzip binary
        bcftools_folder : str, Optional
                          Path to folder containing the bcftools binary
        bedtools_folder : str, Optional
                          Path to the folder containing the bedtools binary
        picard_folder : str, Optional
                        Path to folder containing the picard binary
        r_folder : str, Optional
                   Path to folder containing the R binary
        r_scripts : str, Optional
                    Path to folder containing the R scripts required for
                    constructing some of the plots used by this class (i.e. plot_variants.R)
        tabix_folder : str, Optional
                        Path to folder containing the tabix binary
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
                         outdir=None, intervals=None):
        '''
        Method to calculate the genotype concordance between VcfQC.vcf and Truth VCF.
        It will run Picard's GenotypeConcordance

        Parameters
        ----------
        truth_vcf : str, Required
                    The VCF containing the truth sample
        truth_sample : str, Required
                       The name of the truth sample within the truth VCF
        call_sample : str, Required
                      The name of the call sample within the call VCF
        outprefix :str, Required
                   String used as the prefix in the output file
        outdir : str, optional
                 If provided, then put output files in this folder
        intervals : str, Required
                    One or more interval list files that will be used to limit the
                    genotype concordance

        Returns
        -------
        It returns a GTPconcordance object

        '''

        if outdir:
            outprefix = "%s/%s" % (outdir, outprefix)

        command = "java -jar "
        if self.picard_folder:
            command += self.picard_folder+"/"

        command += "picard.jar GenotypeConcordance TRUTH_VCF={0} CALL_VCF={1} "\
        "TRUTH_SAMPLE={2} CALL_SAMPLE={3} O={4}".format(truth_vcf, self.vcf,
                                                        truth_sample, call_sample, outprefix)

        if intervals:
            command += " INTERVALS={0}".format(intervals)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            raise Exception(exc.output)

        gtp_con = GTPconcordance(summary_metrics_file=outprefix+\
                                 ".genotype_concordance_summary_metrics")

        return gtp_con

    def get_chros(self, filter_str=None, chr_f=None, verbose=None):
        '''
        Method to get a list of chromosomes present in a file

        Parameters
        ----------
        filter_str : str, Optional
                     If defined, apply this filter string so bcftools view
                     apply it before fetching the chros.
        chr_f : str, Optional
                      Path to file with a list of chromosomes (one per line).
                      If provided, the chros in the file will be compared with the
                      chromosomes in self.vcf
        verbose : boolean, optional

        Returns
        -------
        A dict with a key named 'in_vcf' and whose values are the chros that are present in self.vcf

        If list_of_chros is defined, then it will also add 3 keys to the dict:
            'both' whose values will be the chros present in self.vcf and in 'chr_f'
            'in_A' whose values will be the chros PRESENT in self.vcf and NOT in 'chr_f'
            'in_B' whose values will be the chros NOT present in self.vcf and PRESENT in 'chr_f'
        '''

        command = ""
        if self.bcftools_folder:
            command += self.bcftools_folder+"/"

        command += "bcftools view --no-header {0} ".format(self.vcf)

        if filter_str != None:
            command += "-f {0} ".format(filter_str)

        command += "|cut -f1 |uniq"

        out_str = ""

        try:
            if verbose is True:
                print("Command is: %s" % command)
            out = subprocess.check_output(command, shell=True)
            out_str = out.decode("utf-8")
            out_str = out_str.rstrip('\n')

        except subprocess.CalledProcessError as exc:
            raise Exception(exc.output)

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

    def number_variants_in_region(self, region, outprefix):
        '''
        Method to get the number of variants in a particular region/s

        Parameters
        ----------
        region : string, Required
                 String with path to the BED file containing the regions for
                 which the number will be calculated
        outprefix : str, Required
                   Prefix for outfile

        Returns
        -------
        File with the number of variants for each particular region
        '''

        command = ""

        if self.bedtools_folder:
            command += self.bedtools_folder+"/"

        outprefix = outprefix+".counts"

        command += "bedtools coverage -a {0} -b {1} -counts > {2}".\
        format(region, self.vcf, outprefix)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            raise Exception(exc.output)

        return outprefix

    def run_CollectVariantCallingMetrics(self, outprefix, truth_vcf, intervals=None):

        '''
        Method to run Picard's CollectVariantCallingMetrics on a VcfQC object.

        Parameters
        ----------
        outprefix : str, Required
                   Prefix for outfiles: prefix.variant_calling_detail_metrics
                   and prefix.variant_calling_summary_metrics
        truth_vcf : str, Required
                   Reference VCF file
        intervals : str, Optional
                    Target intervals to restrict analysis to

        Returns
        -------
        This function returns a CollectVCallingMetrics object
        '''

        command = "java -jar "
        if self.picard_folder:
            command += self.picard_folder+"/"

        command += "picard.jar CollectVariantCallingMetrics I={0} O={1} "\
        "DBSNP={2}".format(self.vcf, outprefix, truth_vcf)

        if intervals:
            command += " TI={0}".format(intervals)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            raise Exception(exc.output)

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
        length : int, Required
                      Length in bp of the genomic windows for which the number
                      of variants will be calculated
        genome : str, Required
                      File with genome sizes
        outprefix : str, Required
                      Prefix for output file (png file with variant density plot)
        plot_params : dict, Required
                      Dictionary containing the graphical parameters used for plot_variants.R:
                      {
                       'height': int,
                       'chromName_cex': double
                       }

        Returns
        -------
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
        outpath : str, required
              output path
        filter_str : str, optional. Example:  PASS,.
                  Apply filters when calculating the stats
        region : str, optional. Example: chr20
                 Region used for calculating the stats
        region_file : BED file, optional
                      BED file with the regions that will be analyzed
        verbose : boolean, optional

        Returns
        -------
        A BcftoolsStats object
        '''

        command = ""
        if self.bcftools_folder:
            command += self.bcftools_folder+"/"

        outpath = outpath+".stats"

        command += "bcftools stats {0} ".format(self.vcf)

        if region != None:
            command += "-r {0} ".format(region)

        if region_file != None:
            command += "-R {0} ".format(region_file)

        if filter_str != None:
            command += "-f {0} ".format(filter_str)

        command += "> {0}".format(outpath)

        try:
            if verbose is True:
                print("Command is: %s" % command)
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            raise Exception(exc.output)

        #parse output file
        if os.path.isfile(outpath) is False:
            raise Exception("File with bcftools stats %s does not exist" % outpath)

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

        Class variables
        ---------------
        filename : str, Required
                   Filename of the VCF that was used to run bcftools stats
        summary_numbers : dict, Required
                          Dictionary containing the basic stats. i.e.:
                              number of samples:      1
                              number of records:      1867316
                              .....
        ts_tv : float, Required
                ts/tv ratio
        ts_tv_1stalt : float, Required
                ts/tv (1st ALT)
        no_singleton_snps : int, Required
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

        Class variables
        ---------------
        summary_metrics_file : str, Required
                               Filepath to *.genotype_concordance_summary_metrics generated by
                               Picard's GenotypeConcordance
        summary_metrics_snps : dict, Optional
                               Dict with the summary metrics found in the
                               *.genotype_concordance_summary_metrics file (only for SNPs)
        summary_metrics_indels : dict, Optional
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

        Class variables
        ---------------
        vc_detail_metrics_file : str, Required
                                Filepath to *.variant_calling_detail_metrics
                                generated by Picard's CollectVariantCallingMetrics
        vc_summary_metrics_file : str, Required
                                Filepath to *.variant_calling_summary_metrics
                                generated by Picard's CollectVariantCallingMetrics
        vc_detail_metrics : dict, Optional
                            Dictionary with information in file
                            *.variant_calling_detail_metrics
        vc_summary_metrics : dict, Optional
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
