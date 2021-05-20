'''
Created on 27 Feb 2017

@author: ernesto
'''

import os
import json
import glob
import ast
from collections import namedtuple
from Utils.RunProgram import RunProgram

class GATK(object):
    """
    Class to filter a VCF file using different GATK components
    """

    def __init__(self, vcf, reference, caller='UG', bgzip_folder=None, tabix_folder=None,
                 gatk_folder=None, tmp_dir=None):
        """
        Constructor

        Parameters
        ----------
        vcf : str
              Path to gzipped vcf file.
        caller : str, default='UG'
                 Caller used to generate the VCF: UG, bcftools.
        reference : str
                    Path to fasta file containing the reference.
        bgzip_folder : str, optional
                       Path to folder containing the bgzip binary.
        tabix_folder : str, optional
                       Path to folder containing the tabix binary.
        gatk_folder : str, optional
                      Path to folder containing GATK's jar file.
        tmp_dir : str, optional
                  Path to java temporary directory. This needs to be
                  set for GATK modules such as ApplyRecalibrator that
                  fails because there is not enough space in the default java tmp dir.
        """

        if os.path.isfile(vcf) is False:
            raise Exception("File does not exist")

        self.vcf = vcf
        self.caller = caller
        self.reference = reference
        self.bgzip_folder = bgzip_folder
        self.tabix_folder = tabix_folder
        self.gatk_folder = gatk_folder
        self.tmp_dir = tmp_dir

    def run_variantrecalibrator(self, resources, mode,
                                max_gaussians=None, intervals=None,
                                annotations=None, tranches=None,
                                outprefix="recalibrate", verbose=None,
                                log_file=None):
        """
        Run GATK's VariantRecalibrator on a VcfFilter object

        Parameters
        ----------
        resources : str
                    JSON file with resources to add using the -resources option.
        mode : {'SNP','INDEL'}
               Recalibration mode to employ.
        intervals :  chr1:1-1000, optional
                     One or more genomic intervals over which to operate.
        max_gaussians : int, optional
                        Max number of Gaussians for the positive model.
        annotations : list, optional
                      List of annotations to be used. Default=['DP', 'QD', 'FS', 'SOR',
                      'MQ', 'MQRankSum', 'ReadPosRankSum', 'InbreedingCoeff'].
        tranches : list, default=[100.0,99.9,99.0,90.0]
                   Each element in the list will correspond to the  levels of truth
                   sensitivity at which to slice the data. (in percent, that is 1.0
                   for 1 percent).
        outprefix : str, default='recalibrate'
                    out prefix used for -recalFile, -tranchesFile, -rscriptFile.
        verbose : bool, optional
                  Increase verbosity
        log_file : str, optional
                   Path to file that will used for logging the GATK stderr and stdout.

        Returns
        -------
        dict
            Dictionary with location of tranches and recal files
        """

        if annotations is None:
            annotations = ['DP', 'QD', 'FS', 'SOR', 'MQ', 'MQRankSum',
                           'ReadPosRankSum', 'InbreedingCoeff']

        if tranches is None:
            tranches = [100.0, 99.9, 99.0, 90.0]

        if self.caller != 'UG':
            raise Exception("VCF caller %s is incompatible" % self.caller)

        if mode not in ('SNP', 'INDEL'):
            raise Exception("VariantRecalibrator mode is not valid."
                            "Valid values are 'SNP','INDEL'" % mode)

        # prepare the prefix used for output files
        outprefix += "_%s" % mode

        Arg = namedtuple('Argument', 'option value')

        args = [Arg('-T', 'VariantRecalibrator'), Arg('-R', self.reference),
                Arg('-input', self.vcf), Arg('-mode', mode),
                Arg('-recalFile', "{0}.recal".format(outprefix)),
                Arg('-tranchesFile', "{0}.tranches".format(outprefix)),
                Arg('-rscriptFile', "{0}_plots.R".format(outprefix))]

        # Read-in the different resources from the resource JSON file
        resources_str = ""

        with open(resources) as data_file:
            data = json.load(data_file)
            bits = data['resources']
            for dummy, dic in enumerate(bits):
                args.extend([Arg("-resource:%s,known=%s,training=%s,"
                                 "truth=%s,prior=%.1f" % (dic['resource'],
                                                          str(dic['known']).lower(),
                                                          str(dic['training']).lower(),
                                                          str(dic['truth']).lower(),
                                                          dic['prior']), dic['path'])])

        # prepare the -an options
        for elt in annotations:
            args.append(Arg('-an', elt))

        # prepare the list of -tranche option
        if type(tranches) == str:
            tranches = ast.literal_eval(tranches)

        for elt in tranches:
            args.append(Arg('-tranche', elt))

        runner = RunProgram(program='java -jar {0}/GenomeAnalysisTK.jar'.
                            format(self.gatk_folder), args=args,
                            log_file=log_file)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout, stderr, is_error = runner.run_popen()

        recal_f = glob.glob("{0}*.recal".format(outprefix))
        tranches_f = glob.glob("{0}*.tranches".format(outprefix))

        if not recal_f:
            raise Exception("No *.recal files were retrieved after running VariantRecalibrator")
        elif not tranches_f:
            raise Exception("No *.tranches files were retrieved after running VariantRecalibrator")

        if len(recal_f) > 1:
            raise Exception("More than one *.recal file was retrieved")
        elif len(tranches_f) > 1:
            raise Exception("More than one *.tranches file was retrieved")

        return {
            'recal_f': recal_f[0],
            'tranches_f': tranches_f[0]
            }

    def run_applyrecalibration(self, mode, recal_file, tranches_file, outprefix,
                               ts_filter_level=99.0, num_threads=1, compress=True,
                               verbose=None, log_file=None):
        """
        Run GATK's ApplyRecalibration on a VcfFilter object

        Parameters
        ----------
        mode : {'SNP','INDEL'}
               Recalibration mode to employ.
        recal_file : str
                     The input recal file used by ApplyRecalibration.
        tranches_file : str
                        The input tranches file describing where to cut the data.
        outprefix : str
                    Prefix used for the output.
        ts_filter_level : float, default=99.0
                          The truth sensitivity level at which to start filtering.
        num_threads : int, default=1
                      Number of data threads to allocate to this analysis.
        compress : bool, default=True
                   Compress the recalibrated VCF.
        verbose : bool, optional
                  Increase verbosity.
        log_file : str, optional
                   Path to file that will used for logging the GATK stderr and stdout.

        Returns
        -------
        outfile : str
                 Path to filtered VCF file.
        """

        if self.caller != 'UG':
            raise Exception("VCF type %s is incompatible" % self.caller)

        if mode != 'SNP' and mode != 'INDEL':
            raise Exception("ApplyRecalibration mode is not valid."
                            "Valid values are 'SNP','INDEL'" % mode)

        # generate output file name
        outfile = ""
        if mode == 'SNP':
            outfile += "%s.recalibrated_snps_raw_indels.vcf" % outprefix
        elif mode == 'INDEL':
            outfile += "%s.recalibrated_variants.vcf" % outprefix

        Arg = namedtuple('Argument', 'option value')

        args = []

        args.extend([Arg('-jar', '{0}/GenomeAnalysisTK.jar'.format(self.gatk_folder)),
                     Arg('-T', 'ApplyRecalibration'),
                     Arg('-R', self.reference), Arg('-input', self.vcf), Arg('-mode', mode),
                     Arg('--ts_filter_level', ts_filter_level), Arg('-recalFile', recal_file),
                     Arg('--num_threads', num_threads),
                     Arg('-tranchesFile', tranches_file)])

        pipelist = None
        if compress is True:
            outfile += ".gz"
            compressRunner = RunProgram(path=self.bgzip_folder, program='bgzip',
                                        parameters=['-c', '>', outfile])
            pipelist = [compressRunner]
        else:
            args.append(Arg('-o', outfile))

        program_str = None
        if self.tmp_dir is not None:
            program_str = "java -Djava.io.tmpdir={0}".format(self.tmp_dir)
        else:
            program_str = "java"

        runner = RunProgram(program=program_str, args=args, downpipe=pipelist, log_file=log_file)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout, stderr, is_error = runner.run_popen()

        # create an index for the recalibrated file
        if compress is True:
            tabixRunner = RunProgram(path=self.tabix_folder, program='tabix', parameters=[outfile])
            stdout = tabixRunner.run_checkoutput()

        return outfile
