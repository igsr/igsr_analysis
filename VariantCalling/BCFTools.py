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
    '''

    def __init__(self, bam, reference, settings):
        '''
        Constructor

        Class variables
        ---------------
        bam : str, Required
             Path to BAM file used for the variant calling process
        reference : str, Required
             Path to Fasta file containing the reference
        settings : str, Required
                   Path to .ini file with settings
        '''

        if os.path.isfile(bam) is False:
            raise Exception("BAM file does not exist")

        self.bam = bam
        self.reference = reference

        # parse settings file (in .ini file)
        parser = ConfigParser(strict=False)
        parser.optionxform = str

        parser.read(settings)
        self.settings = parser

    def run_bcftools(self, outprefix, annots=['DP', 'SP', 'AD'], m_pileup=1, m_call=False,
                     v=False, O='z', ploidy="GRCh38", threads=1,
                     S=None, r=None, verbose=True):
        '''
        Run BCFTools mpileup and then pipe to BCTools call in order to do the variant calling

        Parameters
        ----------

        outprefix : str, Required
                    Prefix for output VCF file. i.e. /path/to/file/test
        annots : list, Optional
                 mpileup parameter
                 Comma separated list of annotations used to decorate the VCF
        m_mpileup : int, Optional
                    mpileup parameter
                    Minimum number gapped reads for indel candidates. Default=1
        m_call : boolean, Optional
                 call parameter
                 alternative modelfor multiallelic and rare-variant calling
                 designed to overcome known limitations in -c calling model
        v : bool, Optional
            call parameter
            output variant sites only
        O : str, Optional
            call parameter
            output type. Default= 'z'
            Possible values are: BCF (b), uncompressed BCF (u), compressed VCF (z),
            uncompressed VCF (v)
        ploidy : str, Optional
                 predefined ploidy. Default: GRCh38
        threads : int, Optional
                  Number of extra output compression threads.Default=1
        S : str, Optional
            call parameter
            File of sample names to include or exclude if prefixed with "^"
        r: str, Optional
           Region used for doing the variant calling in the format chr20:10000-20000

        verbose : bool, Optional
                  Increase verbosity. Default= True

        Returns
        -------
        A VCF file with variants

        '''

        Arg = namedtuple('Argument', 'option value')

        arguments_mpileup = [Arg('-f', self.reference)]

        settings_mpileup = self.settings.items("bcftools")

        pdb.set_trace()
        for key, v in settings_mpileup:
            arguments_mpileup.append(Arg(key, v))

        arguments_mpileup.append(Arg('-d', d))
        arguments_mpileup.append(Arg('-m', m_pileup))
        arguments_mpileup.append(Arg('--threads', threads))

        for a in annots:
            arguments_mpileup.append(Arg('-a', a))


        if r is not None:
            region_str = re.sub(':|-', '_', r)
            outprefix += ".{0}".format(region_str)
            arguments_mpileup.append(Arg('-r', r))

        pdb.set_trace()
        params_mpileup = []
        if self.settings.has_option('bcftools', 'E'):
            if self.settings.getboolean('bcftools', 'E') is True:
                params_mpileup.append('-E')

        if self.settings.has_option('bcftools', 'p'):
            if self.settings.getboolean('bcftools', 'E') is True:
                params_mpileup.append('-E')
        if p is True:
            params_mpileup.append('-p')

        params_mpileup.append(self.bam)
        params_call = []
        if m_call is True:
            params_call.append('-m')
        if v is True:
            params_call.append('-v')

        arguments_call = []
        arguments_call.append(Arg('-O', O))
        arguments_call.append(Arg('--ploidy', ploidy))
        if S is not None:
            arguments_call.append(Arg('-S', S))

        outprefix += ".vcf.gz"
        arguments_call.append(Arg('-o', outprefix))

        pipelist = None
        bcftools_callRunner = RunProgram(program='bcftools call',
                                         args=arguments_call,
                                         parameters=params_call)
        pipelist = [bcftools_callRunner]

        runner = RunProgram(program='bcftools mpileup',
                            args=arguments_mpileup,
                            parameters=params_mpileup,
                            downpipe=pipelist)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout, stderr, is_error = runner.run_popen(raise_exc=False)

        return outprefix
