'''
Created on 27 Feb 2017

@author: ernesto
'''

import os
import json
import glob
import ast
import pdb
from collections import namedtuple
from typing import List
from Utils.RunProgram import RunProgram

class GATK(object):
    """
    Class to filter a VCF file using different GATK components

    Class variables
    ---------------
    gatk_folder : str, Optional
                  Path to folder containing the gatk3 binary
    bgzip_folder : str, Optional
                   Path to folder containing the bgzip binary
    tabix_folder : str, Optional
                   Path to folder containing the bgzip binary
    arg : namedtuple
          Containing a particular argument and its value
    """
    gatk_folder = None
    bgzip_folder = None
    tabix_folder = None
    arg = namedtuple('Argument', 'option value')

    def __init__(self, vcf : str, caller : str='UG')->None:
        """
        Constructor

        Parameters
        ----------
        vcf : Path to gzipped vcf file.
        caller : Caller used to generate the VCF: UG, bcftools.
        """

        if os.path.isfile(vcf) is False:
            raise Exception("File does not exist")

        self.vcf = vcf
        self.caller = caller

    def run_variantrecalibrator(self, resources: str,
                                annotations: List[str]=None, tranches: List[float]=None,
                                outprefix: str="recalibrate", verbose: bool=False,
                                log_file: str=None, **kwargs)->dict:
        """
        Run GATK's VariantRecalibrator on a VcfFilter object

        Parameters
        ----------
        resources : JSON file with resources to add using the -resources option.
        annotations : Aannotations to be used. Default=['DP', 'QD', 'FS', 'SOR',
                      'MQ', 'MQRankSum', 'ReadPosRankSum', 'InbreedingCoeff'].
        tranches : Each element in the list will correspond to the  levels of truth
                   sensitivity at which to slice the data. (in percent, that is 1.0
                   for 1 percent).
        outprefix : out prefix used for -recalFile, -tranchesFile, -rscriptFile.
        verbose : Increase verbosity.
        log_file : Path to file that will used for logging the GATK stderr and stdout.
        **kwargs: Arbitrary keyword arguments. Check the `GATK` help for
                  more information.

        Returns
        -------
        Dictionary with location of tranches and recal files
        """

        if annotations is None:
            annotations = ['DP', 'QD', 'FS', 'SOR', 'MQ', 'MQRankSum',
                           'ReadPosRankSum', 'InbreedingCoeff']

        if tranches is None:
            tranches = [100.0, 99.9, 99.0, 90.0]

        if self.caller != 'UG':
            raise Exception("VCF caller %s is incompatible" % self.caller)

        if kwargs['mode']:
            # prepare the prefix used for output files
            outprefix += f"_{kwargs['mode']}"

        args = [GATK.arg('-input', self.vcf), 
                GATK.arg('-recalFile', f"{outprefix}.recal"),
                GATK.arg('-tranchesFile', f"{outprefix}.tranches"),
                GATK.arg('-rscriptFile', f"{outprefix}_plots.R")]

        # add now the **kwargs        
        allowed_keys = ['R', 'mode' ] # allowed arbitrary args
        args.extend([GATK.arg(f"-{k}", v) for k, v in kwargs.items() if k in allowed_keys])

        # Read-in the different resources from the resource JSON file
        resources_str = None

        with open(resources) as data_file:
            data = json.load(data_file)
            bits = data['resources']
            for dummy, dic in enumerate(bits):
                args.extend([GATK.arg("-resource:%s,known=%s,training=%s,"
                                      "truth=%s,prior=%.1f" % (dic['resource'],
                                      str(dic['known']).lower(),
                                      str(dic['training']).lower(),
                                      str(dic['truth']).lower(),
                                      dic['prior']), dic['path'])])

        # prepare the -an options
        for elt in annotations:
            args.append(GATK.arg('-an', elt))

        # prepare the list of -tranche option
        if type(tranches) == str:
            tranches = ast.literal_eval(tranches)

        for elt in tranches:
            args.append(GATK.arg('-tranche', elt))

        cmd =  f"{GATK.gatk_folder}/gatk3 -T VariantRecalibrator" if GATK.gatk_folder else "gatk3 -T VariantRecalibrator"
        runner = RunProgram(program=cmd,
                            args=args,
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

    def run_applyrecalibration(self, outprefix: str, log_file: str, compress: bool=True,
                               verbose: bool=True, **kwargs)->str:
        """
        Run GATK's ApplyRecalibration on a VcfFilter object

        Parameters
        ----------
        outprefix : Prefix used for the output.
        compress : Compress the recalibrated VCF.
        log_file : Path to file that will used for logging the GATK stderr and stdout.
        verbose : Increase verbosity.
        **kwargs: Arbitrary keyword arguments. Check the `GATK` help for
                  more information.

        Returns
        -------
        outfile : Path to filtered VCF file.
        """
        if self.caller != 'UG':
            raise Exception("VCF type %s is incompatible" % self.caller)

        args = [GATK.arg('-input', self.vcf)]

         # add now the **kwargs  
        allowed_keys = ['R', 'mode', 'ts_filter_level', 'recalFile', 'nt', 'tranchesFile' ] # allowed arbitrary args
        args.extend([GATK.arg(f"-{k}", v) for k, v in kwargs.items() if k in allowed_keys])

        # generate output file name
        if 'mode' in kwargs:
            if kwargs['mode'] == 'SNP':
                outprefix += f"{outprefix}.recalibrated_snps_raw_indels.vcf"
            elif kwargs['mode'] == 'INDEL':
                outprefix += f"{outprefix}.recalibrated_variants.vcf"

        pipelist = None
        if compress is True:
            outprefix += ".gz"
            cmd_cmp =  f"{GATK.bgzip_folder}/bgzip" if GATK.bgzip_folder else "bgzip"
            compressRunner = RunProgram(program=cmd_cmp,
                                        parameters=['-c', '>', outprefix])
            pipelist = [compressRunner]
        else:
            args.append(GATK.arg('-o', outprefix))

        cmd =  f"{GATK.gatk_folder}/gatk3 -T ApplyRecalibration" if GATK.gatk_folder else "gatk3 -T ApplyRecalibration"
        runner = RunProgram(program=cmd, args=args, downpipe=pipelist, log_file=log_file)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout, stderr, is_error = runner.run_popen()

        # create an index for the recalibrated file
        if compress is True:
            cmd_tabix =  f"{GATK.tabix_folder}/tabix" if GATK.tabix_folder else "tabix"
            tabixRunner = RunProgram(program=cmd_tabix, parameters=[outprefix])
            stdout = tabixRunner.run_checkoutput()

        return outprefix
