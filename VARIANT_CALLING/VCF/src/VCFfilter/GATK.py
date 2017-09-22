'''
Created on 27 Feb 2017

@author: ernesto
'''

import os
import subprocess
import json
import glob
import ast

class GATK(object):
    '''
    Class to filter a VCF file using different GATK components
    '''


    def __init__(self, vcf, caller=None, reference=None, bgzip_folder=None, tabix_folder=None, gatk_folder=None):
        '''
        Constructor

        Class variables
        ---------------
        vcf : str, Required
             Path to gzipped vcf file
        caller : str, Optional
             Caller used to generate the VCF: UG, bcftools
        reference : str, Optional
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
                                outprefix="recalibrate"):
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


        # Read-in the different resources from the resource JSON file
        resources_str = ""

        with open(resources) as data_file:
            data = json.load(data_file)
            bits = data['resources']
            for dummy, dic in enumerate(bits):
                resources_str += "-resource:%s,known=%s,training=%s,truth=%s,prior=%.1f %s " \
                %(dic['resource'], str(dic['known']).lower(), str(dic['training']).lower(), \
                str(dic['truth']).lower(), dic['prior'], dic['path'])

        command = "java -jar "
        if self.gatk_folder:
            command += self.gatk_folder + "/"

        # prepare the -an options
        prefix1 = '-an '
        newlist = [prefix1 + elt for elt in annotations]
        ann_str = " ".join(newlist)

        # prepare the list of -tranche option
        if type(tranches) == str:
            tranches = ast.literal_eval(tranches)
        prefix2 = '-tranche '
        newlist = [prefix2 + str(elt) for elt in tranches]
        tranches_str = " ".join(newlist)

        # prepare the prefix used for output files
        outprefix += "_%s" % mode

        command += "GenomeAnalysisTK.jar -T VariantRecalibrator -R {0} " \
        "-input {1} {2} {3} -mode {4} {5} -recalFile {6}.recal -tranchesFile" \
        " {6}.tranches -rscriptFile {6}_plots.R".format(self.reference, self.vcf, \
                                                         resources_str, ann_str, \
                                                         mode, tranches_str, \
                                                         outprefix)

        if intervals:
            command += " -L {0}".format(intervals)

        if max_gaussians:
            command += " --maxGaussians {0}".format(max_gaussians)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            print("Something went wrong while running VariantRecalibrator.\n"
                  "Command used was: %s" % command)
            raise Exception(exc.output)

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
                               ts_filter_level=99.0, compress=True):
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
        compress : boolean, Default= True
                   Compress the recalibrated VCF
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

        command = "java -jar "
        if self.gatk_folder:
            command += self.gatk_folder + "/"

        if compress is True:
            outfile += ".gz"
            command += "GenomeAnalysisTK.jar -T ApplyRecalibration -R {0} -input {1} -mode {2} \
            --ts_filter_level {3} -recalFile {4} -tranchesFile {5} | {6}/bgzip -c > {7}"\
            .format(self.reference, self.vcf, mode, ts_filter_level, recal_file, tranches_file,\
                     self.bgzip_folder, outfile)
        else:
            command += "GenomeAnalysisTK.jar -T ApplyRecalibration -R {0} -input {1} -mode {2} \
            --ts_filter_level {3} -recalFile {4} -tranchesFile {5} -o {6}"\
            .format(self.reference, self.vcf, mode, ts_filter_level, recal_file, tranches_file,
                    outfile)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exp:
            print("Something went wrong while running ApplyRecalibration")
            raise Exception(exp.output)

        # create an index for the recalibrated file
        if compress is True:
            command = "{0}/tabix {1}".format(self.tabix_folder, outfile)
            try:
                subprocess.check_output(command, shell=True)
            except subprocess.CalledProcessError as exc:
                print("Something went wrong while trying to create the index for the \
                recalibrated file")
                raise Exception(exc.output)

        return outfile
