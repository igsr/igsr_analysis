'''
Created on 16 Aug 2018

@author: ernesto
'''

import os
import re
import pdb
from Utils.RunProgram import RunProgram
from collections import namedtuple

class BCFTools(object):
    """
    Class to run the BCFTools variant caller

    Class variables
    ---------------
    bcftools_folder : str, Optional
                      Path to folder containing the bcftools binary
    arg : namedtuple
          Containing a particular argument and its value
    """
    bcftools_folder = None 
    arg = namedtuple('Argument', 'option value')

    def __init__(self, bam: str, reference: str)->None:
        """
        Constructor

        Parameters
        ----------
        bam : Path to BAM file used for the variant calling process.
        reference : Path to fasta file containing the reference.
        """
        if os.path.isfile(bam) is False:
            raise Exception(f"{bam} file does not exist")

        self.bam = bam
        self.reference = reference

    def run_bcftools(self, *args, prefix: str, annots=['DP', 'SP', 'AD'], threads: int=1,
                     r: str=None, ploidy: str='GRCh38', verbose: bool= False, **kwargs )->str:
        """
        Run BCFTools mpileup and then pipe to BCTools call in order to do the variant calling

        Parameters
        ----------
        prefix : Prefix for output VCF file.
        threads : Number of extra output compression threads.
        r: Region used for doing the variant calling in the format chr20:10000-20000.
        ploidy: Predefined ploidy.
        verbose : Increase verbosity.
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments. Check the `bcftools` help for
                  more information.

        Returns
        -------
        outprefix : VCF file with variants.
        """
        ## mpileup
        arguments_mpileup = [BCFTools.arg('-f', self.reference), BCFTools.arg('--threads', threads) ]

        for a in annots:
            arguments_mpileup.append(BCFTools.arg('-a', a))

        # add now the **kwargs        
        allowed_keys = ['P', 'F', 'C', 'd', 'm'] # allowed arbitrary args
        arguments_mpileup.extend([BCFTools.arg(f"-{k}", v) for k, v in kwargs.items() if k in allowed_keys])

        if r is not None:
            region_str = re.sub(':|-', '_', r)
            prefix += ".{0}".format(region_str)
            arguments_mpileup.append(BCFTools.arg('-r', r))
    
        # add now the *args
        params_mpileup = [self.bam]
        allowed_mpileup_opts = ['E', 'p'] # allowed options
        params_mpileup.extend([f"-{x}" for x in args if x in allowed_mpileup_opts])

        ## call
        # add now the **kwargs 
        allowed_keys = ['O', 'S'] # allowed arbitrary args
        arguments_call = [BCFTools.arg('--ploidy', ploidy)]
        arguments_call.extend([BCFTools.arg(f"-{k}", v) for k, v in kwargs.items() if k in allowed_keys])

        # add now the *args
        allowed_call_opts = ['m', 'v'] # allowed options
        params_call = [f"-{x}" for x in args if x in allowed_call_opts]

        prefix += ".vcf.gz"
        arguments_call.append(BCFTools.arg('-o', prefix))

        pipelist = None
        cmd1 = f"{BCFTools.bcftools_folder}/bcftools call" if BCFTools.bcftools_folder else "bcftools call"
        bcftools_callRunner = RunProgram(program=cmd1,
                                         args=arguments_call,
                                         parameters=params_call)
        pipelist = [bcftools_callRunner]

        cmd2 =  f"{BCFTools.bcftools_folder}/bcftools mpileup" if BCFTools.bcftools_folder else "bcftools mpileup"
        runner = RunProgram(program=cmd2,
                            args=arguments_mpileup,
                            parameters=params_mpileup,
                            downpipe=pipelist)
        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout, stderr, is_error = runner.run_popen(raise_exc=False)

        return prefix
