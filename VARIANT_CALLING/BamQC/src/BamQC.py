'''
Created on 12 Oct 2016

@author: ernesto
'''
from types import *
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import re
import os
import sys
import numpy as np
import pandas as pd
import subprocess
from collections import defaultdict, OrderedDict
import contextlib

import pdb

class BamQC(object):
    '''
    Class to do quality assessment on a BAM file
    '''

    def __init__(self, bam, samtools_folder=None, java_folder=None,
                 picard_folder=None, chk_indel_folder=None,
                 verifybamid_folder=None):
        '''
        Constructor

         Class variables
        ---------------
        bam : str
            Path to BAM file

        samtools_folder : str
            Path to folder containing the samtools binary
        java_folder : str, optional
            Path to folder containing the java binary
        picard_folder : str, optional
            Path to folder containing the Picard jar file
        chk_indel_folder : str, optional
            Path to folder containing Heng Li's chk_indel_rg binary
        verifybamid_folder : str, optional
            Path to folder containing VerifyBAMID
        '''

        if os.path.isfile(bam) is False:
            raise Exception("File does not exist")
        self.bam = bam
        self.java_folder = java_folder
        self.samtools_folder = samtools_folder
        self.picard_folder = picard_folder
        self.chk_indel_folder = chk_indel_folder
        self.verifybamid_folder = verifybamid_folder

    def get_contigs(self):
        '''
        Get all contigs from this BAM

        Parameters
        ----------
        None

        Returns
        -------
        A dictionary containing the following information:
        
            {'contig Name': length (in bp)}
        '''
        header, err = subprocess.Popen(["samtools", "view", "-H", self.bam],
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       env=dict(os.environ, PATH="%s" % self.samtools_folder))\
                                       .communicate()

        header = header.decode("utf-8")
        err = err.decode("utf-8")

        if err != "":
            raise Exception(err)
        # Extract contigs from header and convert contigs to integers
        contigs = {}
        for x in re.findall("@SQ\\WSN:(?P<chrom>[A-Za-z0-9_]*)\\WLN:(?P<length>[0-9]+)", header):
            contigs[x[0]] = int(x[1])
        return contigs

    def list_of_samples(self):
        '''
        Get the samples names from the header of the BAM file

        Parameters
        ----------
        None

        Returns
        -------
        List with the sample names
        '''
        header, err = subprocess.Popen(["samtools", "view", "-H", self.bam], stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       env=dict(os.environ, PATH="%s" % self.samtools_folder))\
                                       .communicate()

        header = header.decode("utf-8")
        err = err.decode("utf-8")

        if err != "":
            raise Exception(err)
        
        samples = re.findall("SM:([\w.]+)\s", header)
        samples = list(set(samples))
        return samples

    def list_of_readgroups(self):
        '''
        Get the Read Groups extracted from the header of the BAM file

        Parameters
        ----------
        None

        Returns
        -------
        List composed of the read groups
        '''
        readgroups = []

        header, err = subprocess.Popen(["samtools", "view", "-H", self.bam],
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       env=dict(os.environ,
                                                PATH="%s" % self.samtools_folder)).communicate()

        header = header.decode("utf-8")
        err = err.decode("utf-8")

        if err != "":
            raise Exception(err)
        for e in header.split('\n'):
            if e.startswith("@RG"):
                readgroups.append(re.findall("ID:([\w.]+)\s", e)[0])

        return list(set(readgroups))

    def get_simple_stats(self):
        '''
        Get a dict with stats on the BAM file as calculated by samtools flagstat

        Parameters
        ----------
        None

        Returns
        -------
        A dictionary containing the following information:
             {
             "total_no_reads": int
             "no_duplicates": int
             "total_no_mapped": int
             "no_properly_paired":  int
             }
        '''
        stats, err = subprocess.Popen(["samtools", "flagstat", self.bam],
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE,
                                      env=dict(os.environ,
                                               PATH="%s" % \
                                               self.samtools_folder)).\
                                               communicate()

        stats = stats.decode("utf-8")
        err = err.decode("utf-8")

        if err != "":
            raise Exception(err)
        stats_list = stats.split("\n")

        stats = {
            "total_no_reads": int(re.search("\\d+", stats_list[0]).group()),
            "no_duplicates": int(re.search("\\d+", stats_list[3]).group()),
            "total_no_mapped": int(re.search("\\d+", stats_list[4]).group()),
            "no_properly_paired":  int(re.search("\\d+", stats_list[8]).group())
            }

        return stats

    #recursive nested dict
    def __rec_dd(self):
        return defaultdict(self.__rec_dd)

    def run_samtools_depth(self, chros):
        '''
        Calculate several coverage metrics on a whole genome sequencing BAM
        file using 'samtools depth'

        Parameters
        ----------
        chros : list or string
            List of contigs or just a single contig used for calculating the coverage

        Returns
        ------
        List of SDepth objects

        This method runs samtools depth on a BAM file and will calculate the following metrics:
            * Number of Bases mapped: This is the number of bases having at least one read mapped
            * Sum of depths of coverage: This is the sum of all the depths in each of the Bases mapped
            * Breadth of coverage: This is the result of dividing bases_mapped/length(contig)
              (i.e. what portion of the contig has reads mapped)
            * Depth of coverage: This is the result of dividing sum_of_depths/length(contig)
        '''

        # Check to see if file exists
        if os.path.isfile(self.bam) is False:
            raise Exception("Bam file does not exist")

        contigs = self.get_contigs()

        if type(chros) is not list:
            chros = [chros]

        list_of_cvgs = []

        for c in chros:
            command = "%s/samtools depth -r %s %s | awk 'BEGIN {max = 0}"\
            "{if ($3>max) max=$3;sum+=$3;cnt++}END{print cnt \"\t\" sum \"\t\" max}'" \
            % (self.samtools_folder, c, self.bam)
            
            bases_mapped, sum_of_depths, max = map(int, subprocess.Popen(command,
                                                                         stdout=subprocess.PIPE,
                                                                         shell=True).\
                                                   communicate()[0].decode("utf-8").strip().split("\t"))
            
            #create Coverage object to hold coverage info
            covO = SDepth()
            covO.contig = c
            covO.max = max
            covO.bases_mapped = bases_mapped
            covO.sum_of_depths = sum_of_depths
            covO.breadth = bases_mapped/float(contigs[c])
            covO.depth = sum_of_depths/float(contigs[c])
            covO.length = int(contigs[c])

            list_of_cvgs.append(covO)

        return list_of_cvgs

    def aggregate_stats(self, cov_list):
        '''
        Used to calculate aggregated stats on a list of SDepth objects

        Parameters
        ----------
        cov_list : list
            List containing the SDepth objects for which the stats will be aggregated.

        Returns
        --------
        A SDepth object

        '''

        assert type(cov_list) is list, "cov_list is not a list: %r" % cov_list

        #create SDepth object to hold coverage info
        covO = SDepth()

        contigstr = ','.join(i.contig for i in cov_list)
        covO.contig = contigstr
        covO.bases_mapped = sum(i.bases_mapped for i in cov_list)
        covO.sum_of_depths = sum(i.sum_of_depths for i in cov_list)
        covO.breadth = sum(i.bases_mapped for i in cov_list)/float(sum(i.length for i in cov_list))
        covO.depth = sum(i.sum_of_depths for i in cov_list)/float(sum(i.length for i in cov_list))
        covO.length = sum(i.length for i in cov_list)
        covO.max = max(i.max for i in cov_list)

        return covO

    def run_CollectHsMetrics(self, baits_file, outfile=None, cov_cap=None):
        '''
        Run Picard's CollectHsMetrics on a Exome sequencing BAM file

        Parameters
        ----------
        baits_file : str, required
            Str consisting on the path to the file containing the Exome baits.
        outfile : str, optional
            If provided, then create a file with the output of this program
        cov_cap : int, optional
            Picard's Coverage Cap parameter. Treat positions with coverage
            exceeding this value as if they had coverage at this value.
            Default value: 250.

        Returns
        ------
        A CMetrics object

        '''

        # Check to see if file exists
        if os.path.isfile(self.bam) is False:
            raise Exception("Bam file does not exist")

        command = ""
        if self.java_folder:
            command += self.java_folder+"/"

        command += "java -jar {0}/picard.jar CollectHsMetrics BI={1} INPUT={2} "\
        "TI={1} OUTPUT=/dev/stdout QUIET=true".format(self.picard_folder, baits_file, self.bam)

        if cov_cap:
            command += " COVERAGE_CAP=%s" % cov_cap

        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stdout, stderr = p.communicate()

        if not stdout:
            raise Exception(stderr)


        if outfile:
            f = open(outfile, 'w')
            f.write(stdout.decode("utf-8"))
            f.close()

        #process stdout
        part = re.split('\n\n', stdout.decode("utf-8"))

        #process metrics
        metrics_keys = part[1].split('\n')[1].split('\t')
        metrics_values = part[1].split('\n')[2].split('\t')

        d = dict(zip(metrics_keys, metrics_values))
        #sort the metrics dictionary by its keys
        sd = OrderedDict(sorted(d.items()))

        #process coverage data to create bar plot
        cov_counts = part[2].split('\n')[1:len(part[2].split('\n'))]
        data = np.array([l.split('\t') for l in part[2].split('\n')[1:len(part[2].split('\n'))]])
        df = pd.DataFrame(data=data[1:, 0:], columns=data[0, 0:])

        df = df.astype(int)
        return CMetrics(metrics=sd, cov_data=df)

    def run_CollectWgsMetrics(self, reference, outfile=None, cov_cap=None):
        '''
        Run Picard's CollectWgsMetrics on a WGS BAM file

        Parameters
        ----------
        reference : str, required
            Str with Fasta file used as the reference.
        outfile : str, optional
            If provided, then create a file with the output of this program
        cov_cap : int, optional
            Picard's Coverage Cap parameter. Treat positions with coverage
            exceeding this value as if they had coverage at this value.
            Default value: 250.

        Returns
        ------
        A CMetrics object

        '''
        if os.path.isfile(self.bam) is False:
            raise Exception("Bam file does not exist")

        command = ""
        if self.java_folder:
            command += self.java_folder+"/"

        command += "java -jar {0}/picard.jar CollectWgsMetrics I={1} OUTPUT="\
        "/dev/stdout R={2} QUIET=true".format(self.picard_folder, self.bam, reference)

        if cov_cap:
            command += " COVERAGE_CAP=%s" % cov_cap

        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stdout, stderr = p.communicate()

        stdout = stdout.decode("utf-8")
        stderr = stderr.decode("utf-8")

        if not stdout:
            raise Exception(stderr)

        if outfile:
            f = open(outfile, 'w')
            f.write(stdout.decode("utf-8"))
            f.close()

        #process stdout
        part = re.split('\n\n', stdout)

        #process metrics
        metrics_keys = part[1].split('\n')[1].split('\t')
        metrics_values = part[1].split('\n')[2].split('\t')

        d = dict(zip(metrics_keys, metrics_values))
        #sort the metrics dictionary by its keys
        sd = OrderedDict(sorted(d.items()))

        #process coverage data to create bar plot
        data = np.array([l.split('\t') for l in part[2].split('\n')[1:len(part[2].split('\n'))]])
        df = pd.DataFrame(data=data[1:, 0:], columns=data[0, 0:])

        df = df.astype(int)
        return CMetrics(metrics=sd, cov_data=df)

    def run_chk_indel_rg(self, outfile=None):
        '''
        Run Heng Li's chk_indel_rg on a BAM file

        Parameters
        ----------
        outfile : str, optional
            If provided, then create a file with the output of this program

        Returns
        ------
        A list of Chk_indel objects

        '''
        if os.path.isfile(self.bam) == False:
            raise Exception("Bam file does not exist")

        command = ""
        if self.chk_indel_folder:
            command += self.chk_indel_folder+"/"

        command += "chk_indel_rg {0}".format(self.bam)

        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stdout, stderr = p.communicate()

        stdout = stdout.decode("utf-8")
        stderr = stderr.decode("utf-8")

        if not stdout:
            raise Exception(stderr)

        if outfile:
            f = open(outfile, 'w')
            f.write(stdout.decode("utf-8"))
            f.close()

        #initialize the list that will contain the Chk_indel objects
        data = []
        stdout = stdout.rstrip('\n')
        for l in stdout.split('\n'):
            data.append(Chk_indel(RG=l.split('\t')[1], ins_in_short_homopolymer=\
                                  int(l.split('\t')[2]), del_in_short\
                                  =int(l.split('\t')[3]), ins_in_long\
                                  =int(l.split('\t')[4]), del_in_long\
                                  =int(l.split('\t')[5])))
        return data

    def run_verifybamid(self, genotype_file, outprefix, outdir=None):
        '''
        Run VerifyBAMID to check for sample swap or contamination issues

        Parameters
        ----------
        genotype_file : str, required
            vcf file with chip genotypes to use
        outprefix : str, required
            prefix for outputfiles
        outdir : str, optional
            If provided, then put output files in this folder

        Returns
        ------
        A list with the paths to the output files generated by VerifyBAMID

        '''

        if os.path.isfile(self.bam) == False:
            raise Exception("Bam file does not exist")

        if outdir: outprefix = "%s/%s" % (outdir, outprefix)

        list_of_outfiles = [outprefix+x for x in ['.depthRG', '.depthSM', '.selfRG', '.selfSM']]

        # check if all output files already exist before running verifybamid
        for s in list_of_outfiles:
            if os.path.isfile(s) == True:
                raise Exception("%s already exists!. It will not overwrite!" % s)

        command = ""
        if self.verifybamid_folder:
            command += self.verifybamid_folder+"/"

        command += "verifyBamID --vcf {0} --bam {1} --out {2} --minAF 0.01"\
        " --minCallRate 0.95 --self".format(genotype_file, self.bam, outprefix)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as e:
            print(e.output)

        # check if all output files exist
        for s in list_of_outfiles:
            if os.path.isfile(s) == False:
                raise Exception("%s was not created!. Check VerifyBAMID run!" % s)

        return list_of_outfiles

class SDepth(object):
    '''
    Class to store coverage metrics on a Whole Genome Sequencing BAM file
    calculated using SAMtools depth
    '''

    def __init__(self, contig=None, mapped=None, breadth=None, depth=None,
                 length=None, sum_of_depths=None, max=None):
        '''
        Create a SDepth object
        '''
        self.contig = contig
        self.bases_mapped = mapped
        self.breadth = breadth
        self.depth = depth
        self.sum_of_depths = sum_of_depths
        self.length = length
        self.max = max

    @contextlib.contextmanager
    def __smart_open(self, filename=None):
        if filename and filename != '-':
            fh = open(filename, 'w')
        else:
            fh = sys.stdout

        try:
            yield fh
        finally:
            if fh is not sys.stdout:
                fh.close()

    def print_report(self, filename=None):
        '''
        Used to print a text report of data in the object

        Parameters
        ----------
        filename : str, optional
            Filename to write the report. The default is STDOUT.
        '''

        with self.__smart_open(filename) as fh:
            print >>fh, 'Stats for %s (length=%d bp):' % (self.contig, self.length)
            print >>fh, '\tSum of depths of coverage: %d' % self.sum_of_depths
            print >>fh, '\tBases mapped: %d' % self.bases_mapped
            print >>fh, '\tDepth of coverage: %.3f' % self.depth
            print >>fh, '\tBreadth of coverage: %.3f' % self.breadth
            print >>fh, '\tMaximum coverage %d' % self.max

class CMetrics(object):
    '''
    Class to store coverage information on the metrics calculated by Picard's
    CollectHsMetrics/CollectWgsMetrics on an Exome or WGS BAM file
    '''

    def __init__(self, metrics, cov_data):
        '''
        Create an CMetrics object

        Class variables
        ---------------
        metrics : dict
            Dictionary with all the metrics generated by running Picard's
            CollectHsMetrics or CollectWgsMetrics
        cov_data : Panda's DataFrame containing data used to generate
                   a bar plot with the coverage counts

        '''
        self.metrics = metrics
        self.cov_data = cov_data

    @contextlib.contextmanager
    def __smart_open(self, filename=None):
        if filename and filename != '-':
            fh = open(filename, 'w')
        else:
            fh = sys.stdout

        try:
            yield fh
        finally:
            if fh is not sys.stdout:
                fh.close()

    def print_report(self, filename=None):
        '''
        Used to print a text report of data in the object

        Parameters
        ----------
        filename : str, optional
            Filename to write the report. The default is STDOUT.
        '''

        with self.__smart_open(filename) as fh:
            keylist = self.metrics.keys()
            keylist.sort()
            for key in keylist:
                print >>fh, key, self.metrics.get(key)

    def create_cov_barplot(self, filename, xlim=None, ylim=None):
        '''
        This method will create a Barplot using the different coverage
        values counts calculated by Picard's
        CollectHsMetrics or CollectWgsMetrics

        Parameters
        ----------
        filename : str, required
            PDF file to write the plot.
            xlim : tuple, optional
            Set the X-axis limit
            ylim : tuple, optional
            Set the Y-axis limit
        '''

        #getting basename from filename
        basename = os.path.basename(filename).split('.')[0]

        ax = None
        if xlim:
            ax = self.cov_data[xlim[0]:xlim[1]].plot(x='coverage', y='count',
                                                     kind="bar", legend=False,
                                                     grid=True, figsize=(20, 10),
                                                     ylim=ylim, color="gold",
                                                     fontsize=14)
        else:
            ax = self.cov_data.plot(x='coverage', y='count', kind="bar",
                                    legend=False, grid=True, figsize=(20, 10),
                                    ylim=ylim, color="gold", fontsize=14)

        plt.title(basename, fontsize=18)
        plt.xlabel('coverage', fontsize=16)
        plt.ylabel('count', fontsize=16)

        ticks = ax.xaxis.get_ticklocs()

        n = 10
        ticklabels = [l.get_text() for l in ax.xaxis.get_ticklabels()]
        ax.xaxis.set_ticks(ticks[::n])
        ax.xaxis.set_ticklabels(ticklabels[::n])
        fig = ax.get_figure()
        fig.savefig(filename, format='pdf')

class Chk_indel(object):
    '''
    Class to store information on the ratio of short insertion and deletion
    calculated by runnint Heng Li's chk_indel_rg
    '''

    def __init__(self, RG, ins_in_short_homopolymer, del_in_short, ins_in_long, del_in_long):
        '''
        Create a Chk_indel object

        Class variables
        ---------------
        RG : readgroup that will be analyssed
        ins_in_short_homopolymer : float, required
                                   ins_in_short_homopolymer
        del_in_short : float, required
                       del_in_short
        ins_in_long : float, required 
                      ins_in_long
        del_in_long : float, required
                      del_in_long

        '''
        self.RG = RG
        self.ins_in_short_homopolymer = ins_in_short_homopolymer
        self.del_in_short = del_in_short
        self.ins_in_long = ins_in_long
        self.del_in_long = del_in_long
