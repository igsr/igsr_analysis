'''
Created on 22 Feb 2018

@author: ernesto
'''

import re
import os
import pdb
from Utils.RunProgram import RunProgram
from collections import namedtuple

class GATK(object):
    """
    Class to run the different GATK variant callers

    Class variables
    ---------------
    gatk_folder : str, Optional
                  Path to folder containing the gatk binary
    bgzip_folder : str, Optional
                   Path to the folder containing the bgzip binary
    arg : namedtuple
          Containing a particular argument and its value
    """
    gatk_folder = None
    bgzip_folder = None
    arg = namedtuple('Argument', 'option value')

    def __init__(self, bam, reference)->None:
        """
        Constructor

        Parameters
        ----------
        bam : str
             Path to BAM file used for the variant calling process.
        reference : str
             Path to fasta file containing the reference.
        """

        if os.path.isfile(bam) is False:
            raise Exception("BAM file does not exist")

        self.bam = bam
        self.reference = reference

    def run_caller(self, program: str, prefix: str, compress: bool=True, log_file: str=None, verbose: bool=False, **kwargs)->str:
        """
        Run GATK UnifiedGenotyper/HaplotypeCaller

        Parameters
        ----------
        program : str
                  Caller used: UnifiedGenotyper, HaplotypeCaller
        prefix : str
                 Prefix for output VCF file. i.e. /path/to/file/test.
        compress : bool, default=True
                   Compress the output VCF.
        log_file : str, optional
                   Path to file that will used for logging the GATK stderr and stdout.
        verbose : bool, default=False
                  Increase verbosity.
        **kwargs: Arbitrary keyword arguments. Check the `GATK` help for
                  more information.

        Returns
        -------
        prefix : str
                A VCF file.
        """

        arguments = [GATK.arg('-T', program), GATK.arg('-R', self.reference),
                     GATK.arg('-I', self.bam)]

        # add now the **kwargs 
        allowed_keys = ['glm', 'nt', 'L', 'alleles', 'gt_mode', 'out_mode'] # allowed arbitrary args
        arguments.extend([GATK.arg(f"-{k}", v) for k, v in kwargs.items() if k in allowed_keys])

        pipelist = None
        if compress is True:
            prefix += ".vcf.gz"
            cmd1 = f"{GATK.bgzip_folder}/bgzip" if GATK.bgzip_folder else "bgzip"
            compressRunner = RunProgram(program=cmd1,
                                        parameters=['-c', '>', prefix])
            pipelist = [compressRunner]
        else:
            prefix += ".vcf"
            arguments.append(GATK.arg('-o', prefix))

        cmd2 = f"{GATK.gatk_folder}/gatk3" if GATK.gatk_folder else "gatk3"
        runner = RunProgram(program=cmd2,
                            args=arguments, downpipe=pipelist, log_file=log_file)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout, stderr, is_error = runner.run_popen(raise_exc=False)

        if is_error is True:
            """
            This piece of code is necessary as GATK crashes when the intersection between
            the genomic chunk and the alleles passed in the VCF
            are calculated and there are no sites.

            If that's the case then GATK  will be run without the interval intersection
            """
            patt = re.compile('##### ERROR MESSAGE: Bad input: '
                              'The INTERSECTION of your -L options produced no intervals.')
            lines = stderr.split('\n')
            interval_error_seen = False
            for l in lines:
                m = patt.match(l)
                if m:
                    interval_error_seen = True
                    alleles = ([arg.value for arg in arguments if arg.option == '--alleles'])[0]
                    for k, i in enumerate(arguments):
                        if i.option == '--intervals' and i.value == alleles:
                            del arguments[k]
                        elif i.option == '--interval_set_rule':
                            del arguments[k]
            if interval_error_seen is False:
                raise Exception(stderr)
            elif interval_error_seen is True:
                runner = RunProgram(program='java -jar {0}/GenomeAnalysisTK.jar'.format(self.gatk_folder),
                                    downpipe=pipelist, args=arguments, log_file=log_file)
                if verbose is True:
                    print("Command line is: {0}".format(runner.cmd_line))

                stdout, stderr, is_error = runner.run_popen(raise_exc=True)

        return prefix