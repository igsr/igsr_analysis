'''
Created on 27 Feb 2017

@author: ernesto
'''

import os
import subprocess
from Utils.RunProgram import RunProgram
from collections import namedtuple
import json
import glob
import ast
import pdb
import re

class GATK(object):
    '''
    Class to filter a VCF file using different GATK components
    '''


    def __init__(self, vcf, reference, caller='UG',  bgzip_folder=None, tabix_folder=None, gatk_folder=None):
        '''
        Constructor

        Class variables
        ---------------
        vcf : str, Required
             Path to gzipped vcf file
        caller : str, Optional
             Caller used to generate the VCF: UG, bcftools
             Default= 'UG'
        reference : str, Required
             Path to fasta file containing the reference
        bgzip_folder : str, Optional
                      Path to folder containing the bgzip binary
        tabix_folder : str, Optional
                        Path to folder containing the tabix binary
        gatk_folder : str, Optional
                       Path to folder containing GATK's jar file
        '''

        if os.path.isfile(vcf) is False:
            raise Exception("File does not exist")

        self.vcf = vcf
        self.caller = caller
        self.reference = reference
        self.bgzip_folder = bgzip_folder
        self.tabix_folder = tabix_folder
        self.gatk_folder = gatk_folder

    def run_variantrecalibrator(self, resources, mode,
                                max_gaussians=None, intervals=None,
                                annotations=None, tranches=None,
                                outprefix="recalibrate",verbose=None):
        '''
        Run GATK's VariantRecalibrator on a VcfFilter object

        Parameters
        ----------
        resources : JSON file with resources to add using the -resources option, Required
        mode : str, Required
                    Recalibration mode to employ (SNP|INDEL)
        intervals :  chr1:1-1000, Optional
                    One or more genomic intervals over which to operate
        max_gaussians : int, Optional
                        Max number of Gaussians for the positive model
        annotations : list, Optional
                      List of annotations to be used. Default=['DP', 'QD', 'FS', 'SOR',
                      'MQ', 'MQRankSum', 'ReadPosRankSum', 'InbreedingCoeff']
        tranches : list, Optional
                   Each element in the list will correspond to the  levels of truth
                   sensitivity at which to slice the data. (in percent, that is 1.0
                   for 1 percent). Default=[100.0,99.9,99.0,90.0]
        outprefix : str, Optional
                    out prefix used for -recalFile, -tranchesFile, -rscriptFile.
                    Default= recalibrate
        verbose : bool, Optional
                  Increase verbosity

        Returns
        -------
        Dictionary with location of tranches and recal files

        '''

        if annotations is None:
            annotations = ['DP', 'QD', 'FS', 'SOR', 'MQ', 'MQRankSum',
                           'ReadPosRankSum', 'InbreedingCoeff']

        if tranches is None:
            tranches = [100.0, 99.9, 99.0, 90.0]

        if self.caller != 'UG':
            raise Exception("VCF caller %s is incompatible" % self.caller)

        if mode != 'SNP' and mode != 'INDEL':
            raise Exception("VariantRecalibrator mode is not valid."
                            "Valid values are 'SNP','INDEL'" % mode)

        # prepare the prefix used for output files
        outprefix += "_%s" % mode

        Arg = namedtuple('Argument', 'option value')

        args=[Arg('-T', 'VariantRecalibrator'), Arg('-R', self.reference), Arg('-input', self.vcf), Arg('-mode', mode),
              Arg('-recalFile', "{0}.recal".format(outprefix)), Arg('-tranchesFile', "{0}.tranches".format(outprefix)),
              Arg('-rscriptFile', "{0}_plots.R".format(outprefix)) ]

        # Read-in the different resources from the resource JSON file
        resources_str = ""

        with open(resources) as data_file:
            data = json.load(data_file)
            bits = data['resources']
            for dummy, dic in enumerate(bits):
                args.extend([Arg("-resource:%s,known=%s,training=%s,truth=%s,prior=%.1f" % (dic['resource'],str(dic['known']).lower(),str(dic['training']).lower(),str(dic['truth']).lower(),dic['prior']), dic['path'])])


        # prepare the -an options
        for elt in annotations:
            args.append(Arg('-an',elt))

        # prepare the list of -tranche option
        if type(tranches) == str:
            tranches = ast.literal_eval(tranches)

        for elt in tranches:
            args.append(Arg('-tranche', elt))

        runner=RunProgram(program='java -jar {0}/GenomeAnalysisTK.jar'.format(self.gatk_folder), args=args)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout,stderr=runner.run_popen()

        lines=stderr.split("\n")
        p = re.compile('#* ERROR')
        for i in lines:
            m = p.match(i)
            if m:
                print("Something went wrong while running GATK VariantRecalibrator. This was the error message: {0}".format(stderr))
                raise Exception()


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
                               ts_filter_level=99.0, num_threads=1, compress=True, verbose=None):
        '''
        Run GATK's ApplyRecalibration on a VcfFilter object

        Parameters
        ----------
        mode : str, Required
                    Recalibration mode to employ (SNP|INDEL)
        recal_file : str, Required
                     The input recal file used by ApplyRecalibration
        tranches_file : str, Required
                        The input tranches file describing where to cut the data
        outprefix : str, Required
                    Prefix used for the output
        ts_filter_level : float, Optional
                          The truth sensitivity level at which to start filtering. Default=99.0
        num_threads : int, Optional
                   Number of data threads to allocate to this analysis. Default=1
        compress : boolean, Default= True
                   Compress the recalibrated VCF
        verbose : bool, Optional
                  Increase verbosity
        '''

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

        args=[Arg('-T', 'ApplyRecalibration'), Arg('-R', self.reference), Arg('-input', self.vcf), Arg('-mode', mode),
              Arg('--ts_filter_level', ts_filter_level), Arg('-recalFile', recal_file), Arg('--num_threads', num_threads),
              Arg('-tranchesFile', tranches_file) ]

        pipelist=None
        if compress is True:
            outfile += ".gz"
            compressRunner=RunProgram(path=self.bgzip_folder,program='bgzip',parameters=[ '-c', '>', outfile])
            pipelist=[compressRunner]
        else:
            args.append(Arg('-o',outfile))

        runner=RunProgram(program='java -jar {0}/GenomeAnalysisTK.jar'.format(self.gatk_folder), args=args, downpipe=pipelist)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout,stderr=runner.run_popen()

        lines=stderr.split("\n")
        p = re.compile('#* ERROR')
        for i in lines:
            m = p.match(i)
            if m:
                print("Something went wrong while running GATK ApplyRecalibration. This was the error message: {0}".format(stderr))
                raise Exception()

        # create an index for the recalibrated file
        if compress is True:
            tabixRunner=RunProgram(path=self.tabix_folder,program='tabix',parameters=[outfile])
            stdout=tabixRunner.run_checkoutput()

        return outfile

