'''
Created on 16 Aug 2018

@author: ernesto
'''

import os
import re
from Utils.RunProgram import RunProgram
from collections import namedtuple

class BCFTools(object):
    """
    Class to run the BCFTools variant caller
    """

    def __init__(self, bam, reference, bcftools_folder):
        """
        Constructor

        Parameters
        ----------
        bam : str
             Path to BAM file used for the variant calling process.
        reference : str
             Path to fasta file containing the reference.
        bcftools_folder : str
                          Path to folder containing the BCFTools binary.
        """

        if os.path.isfile(bam) is False:
            raise Exception("BAM file does not exist")

        self.bam = bam
        self.reference = reference
        self.bcftools_folder = bcftools_folder

    def run_bcftools(self, outprefix, E=False, p=False, annots=['DP', 'SP', 'AD'],
                     P="ILLUMINA", F=0.002, C=50, m_pileup=1, m_call=False,
                     d=250, v=False, O='z', ploidy="GRCh38", threads=1,
                     S=None, r=None, verbose=True):
        """
        Run BCFTools mpileup and then pipe to BCTools call in order to do the variant calling

        Parameters
        ----------
        outprefix : str
                    Prefix for output VCF file. i.e. /path/to/file/test.
        E : bool, default=False
            mpileup parameter
            Recalculate BAQ on the fly, ignore existing BQ tags.
        p : bool, optional
            mpileup parameter
            Apply -m and -F thresholds per sample to increase sensitivity of calling.
            By default both options are applied to reads pooled from all samples.
        annots : list, optional
                 mpileup parameter
                 Comma separated list of annotations used to decorate the VCF.
        P : str, default=ILLUMINA
            mpileup parameter
            Comma-delimited list of patforms (determined by @RG-PL) from which indel
            candidates are obtained.
        F : float, default=0.002
            mpileup parameter
            Minimum fraction of gapped reads.
        C : int, default=50
            mpileup parameter
            Coefficient for downgrading mapping quality for reads containing excessive
            mismatches.
        d : int, default=250
            mpileup parameter
            At a position, read maximally INT reads per input file.
        m_mpileup : int, default=1
                    mpileup parameter
                    Minimum number gapped reads for indel candidates.
        m_call : bool, optional
                 call parameter
                 alternative modelfor multiallelic and rare-variant calling
                 designed to overcome known limitations in -c calling model.
        v : bool, optional
            call parameter
            output variant sites only.
        O : str, default='z'
            call parameter
            output type.
            Possible values are: BCF (b), uncompressed BCF (u), compressed VCF (z),
            uncompressed VCF (v).
        ploidy : str, default=GRCh38
                 predefined ploidy.
        threads : int, default=1
                  Number of extra output compression threads.
        S : str, optional
            call parameter
            File of sample names to include or exclude if prefixed with "^".
        r: str, optional
           Region used for doing the variant calling in the format chr20:10000-20000.
        verbose : bool, default=True
                  Increase verbosity.

        Returns
        -------
        outprefix : str
            A VCF file with variants.
        """

        Arg = namedtuple('Argument', 'option value')


        arguments_mpileup = [Arg('-f', self.reference)]

        for a in annots:
            arguments_mpileup.append(Arg('-a', a))

        arguments_mpileup.append(Arg('-P', P))
        arguments_mpileup.append(Arg('-F', F))
        arguments_mpileup.append(Arg('-C', C))
        arguments_mpileup.append(Arg('-d', d))
        arguments_mpileup.append(Arg('-m', m_pileup))
        arguments_mpileup.append(Arg('--threads', threads))

        if r is not None:
            region_str = re.sub(':|-', '_', r)
            outprefix += ".{0}".format(region_str)
            arguments_mpileup.append(Arg('-r', r))

        params_mpileup = []
        if E is True:
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

        runner = RunProgram(program='{0}/bcftools mpileup'.format(self.bcftools_folder),
                            args=arguments_mpileup,
                            parameters=params_mpileup,
                            downpipe=pipelist)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout, stderr, is_error = runner.run_popen(raise_exc=False)

        return outprefix
