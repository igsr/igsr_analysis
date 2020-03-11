'''
Created on 16 Aug 2018

@author: ernesto
'''

import os
import re
import pdb
from Utils.RunProgram import RunProgram
from collections import namedtuple
from configparser import ConfigParser

class BCFTools(object):
    '''
    Class to run the BCFTools variant caller

    Class Variables
    ---------------
    settings : ConfigParser, Required
               Object with configuration settings
               relevant for this class
    reference : Fasta file containing the reference
                used for alignment
    '''

    def __init__(self, bam, reference_f, settingf):
        '''
        Constructor

        Instance variables
        ------------------
        bam : str, Required
             Path to BAM file used for the variant calling process
        reference_f : str, Required
             Path to Fasta file containing the reference
        settingf : str, Required
                   Path to .ini file with settings
        '''

        if os.path.isfile(bam) is False:
            raise Exception("BAM file does not exist")

        self.bam = bam
        self.reference = reference_f

        # parse settings file (in .ini file)
        parser = ConfigParser(allow_no_value=True)
        parser.optionxform = str

        parser.read(settingf)
        self.settings = parser

    def run_bcftools(self, outprefix, r=None, verbose=True):
        '''
        Run BCFTools mpileup and then pipe to BCTools call in order to do the variant calling

        Parameters
        ----------

        outprefix : str, Required
                    Prefix for output VCF file. i.e. /path/to/file/test

        r: str, Optional
           Region used for doing the variant calling in the format chr20:10000-20000
        verbose : bool, Optional
                  Increase verbosity. Default= True

        Returns
        -------
        A VCF file with variants

        '''

        Arg = namedtuple('Argument', 'option value')

        # mpileup
        mpileup_args_l = [Arg('-f', self.reference)]

        for k, v in self.settings.items("mpileup_opts"):
            mpileup_args_l.append(Arg(k, v))

        if r is not None:
            region_str = re.sub(':|-', '_', r)
            outprefix += ".{0}".format(region_str)
            mpileup_args_l.append(Arg('-r', r))

        mpileup_params_l = []
        for i in self.settings['mpileup_flags'].keys():
            if self.settings['mpileup_flags'][i] is None:
                mpileup_params_l.append(i)

        # add also the BAM file
        mpileup_params_l.append(self.bam)

        # call
        call_args_l = []
        for k, v in self.settings.items("call_opts"):
            call_args_l.append(Arg(k, v))

        outprefix += ".vcf.gz"
        call_args_l.append(Arg('-o', outprefix))

        call_params_l = []
        for i in self.settings['call_flags'].keys():
            if self.settings['call_flags'][i] is None:
                call_params_l.append(i)

        pipelist = None
        bcftools_callRunner = RunProgram(program='bcftools call',
                                         args=call_args_l,
                                         parameters=call_params_l)
        pipelist = [bcftools_callRunner]

        runner = RunProgram(program='bcftools mpileup',
                            args=mpileup_args_l,
                            parameters=mpileup_params_l,
                            downpipe=pipelist)
        pdb.set_trace()
        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout, stderr, is_error = runner.run_popen(raise_exc=False)

        return outprefix
