'''
Created on 21 Jul 2017

@author: ernesto
'''
import os
from collections import namedtuple
from Utils.RunProgram import RunProgram

class Shapeit(object):
    """
    Class to run SHAPEIT
    """

    def __init__(self, bgzip_folder=None, shapeit_folder=None, ligateHAPLOTYPES_folder=None):
        """
        Constructor

        Parameters
        ----------
        bgzip_folder : str, optional
                       Path to folder containing the bgzip binary.
        shapeit_folder : str, optional
                         Path to folder containing the Shapeit binary.
        ligateHAPLOTYPES_folder : str, optional
                                  Path to folder containing the ligateHAPLOTYPES binary.
        """
        self.bgzip_folder = bgzip_folder
        self.shapeit_folder = shapeit_folder
        self.ligateHAPLOTYPES_folder = ligateHAPLOTYPES_folder

    def run_shapeit(self, output_prefix, input_gen=None, input_init=None, input_scaffold=None,
                    input_bed=None, duohmm=False, input_map=None, verbose=False, **kwargs):
        """
        Run Shapeit

        Parameters
        ----------
        input_gen : str, optional
                    Specifies the genotype/GL input data that you obtain from Beagle4,
                    i.e. 'input.shapeit.20.gen.gz input.shapeit.20.gen.sample'.
        input_init : str, optional
                     Specifies the haplotypes that you obtain from Beagle4,
                     i.e. 'input.shapeit.20.hap.gz input.shapeit.20.hap.sample'.
        input_scaffold : str, optional
                         SNP-array derived haplotype scaffold used by SHAPEIT.
                         It has to be in Impute2 format.
                         i.e. 'scaffold.haps.gz scaffold.haps.sample'.
        input_bed : str, optional
                    Unphased genotypes in Plink Binary BED/BIM/FAM format.
                    i.e. 'file.bed file.bim file.fam'.
        duohmm : bool, default=False
                 If true, then activate the --duohmm option.
        output_prefix : str
                        Prefix used for the 2 output files estimated by SHAPEIT,
                        i.e. 'output.shapeit.20.haps.gz output.shapeit.20.haps.sample'.
        input_map : str, optional
                    Path to the file with the genetic map.
        i_from : int, optional
                 Specify the region to be phased.
        i_to : int, optional
               Specify the region to be phased.
        verbose : bool, optional
                  If true, then print the command line used for running this program.

        Returns
        -------
        outdict : dict
            A dict with the path to the 2 output files
            (*.haps.gz and *.haps.sample) that can be used with SHAPEIT.
        """

        if input_gen is None and input_bed is None:
            raise Exception("Error! Either --input-gen or --input-bed need to "
                            "be specified as input for SHAPEIT")

        Arg = namedtuple('Argument', 'option value')

        args = []

        if input_gen is not None:
            args.append(Arg('-call --input-gen', input_gen))
        elif input_bed is not None:
            args.append(Arg('--input-bed', input_bed))

        if input_init is not None:
            args.append(Arg('--input-init', input_init))

        if input_scaffold is not None:
            args.append(Arg('--input-scaffold', input_scaffold))

        if input_map is not None:
            args.append(Arg('--input-map', input_map))

        for k, v in kwargs.items():
            args.append(Arg('--{0}'.format(k), v))

        args.extend([Arg('--output-max',
                         '{0}.haps.gz {0}.haps.sample'.format(output_prefix)),
                     Arg('--output-log', '{0}.log'.format(output_prefix))])

        params = []
        if duohmm is True:
            params = ['--duohmm']

        runner = RunProgram(path=self.shapeit_folder,
                            program='shapeit',
                            args=args,
                            parameters=params)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout = runner.run_checkoutput()

        outdict = {
            'hap_gz': '{0}.haps.gz'.format(output_prefix),
            'hap_sample': '{0}.haps.sample'.format(output_prefix)
        }

        return outdict

    def ligate_shapeitchunks(self, vcf_f, scaffolded_samples, chunk_str,
                             output_prefix, verbose=False):
        """
        Run ligateHAPLOTYPES to ligate together all haplotype chunks produced by SHAPEIT
        (see https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#haplegsample)

        Parameters
        ----------
        vcf_f : str
                VCF file with the Genotype likelihoods.
        scaffolded_samples : str
                             File with the list of samples (separated by '\n')
                             that have been scaffolded.
        chunk_str : str
                    String with the paths to the different files generated by SHAPEIT for
                    the different chromosome chunks
                    (i.e. 's2.chunk1.hap.gz s2.chunk1.hap.gz s2.chunk1.hap.gz').
        output_prefix : str
                        String with the output prefixes (i.e. 'output.shapeit.22.ligated.haps.gz
                        output.shapeit.22.ligated.haps.sample').

        Returns
        -------
        outdict : dict
            A dict with the path to the 2 output files (*.haps.gz and *.haps.sample).
        """

        if self.ligateHAPLOTYPES_folder is None:
            raise Exception("ligateHAPLOTYPES_folder must be defined")

        Arg = namedtuple('Argument', 'option value')

        args = [Arg('--vcf', vcf_f), Arg('--scaffold', scaffolded_samples),
                Arg('--chunks', chunk_str), Arg('--output', '{0}.ligated.haps.gz '
                                                            '{0}.ligated.haps.sample'.
                                                format(output_prefix))]

        runner = RunProgram(path=self.ligateHAPLOTYPES_folder,
                            program='ligateHAPLOTYPES',
                            args=args)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout = runner.run_checkoutput()

        outdict = {
            'hap_gz': '{0}.ligated.haps.gz'.format(output_prefix),
            'hap_sample': '{0}.ligated.haps.sample'.format(output_prefix)
        }

        return outdict

    def convert2vcf(self, input_prefix, output_prefix, compress=False, verbose=False, logfile=None):
        """
        Function to use SHAPEIT's -convert in order to convert the
        *.haps.gz & *.haps.sample files into VCF

        Parameters
        ----------
        input_prefix : str
                       Prefix for the files in HAPS/SAMPLE format.
        output_prefix : str
                        String with the output prefix for the VCF file.
        verbose : bool, optional
                  if true, then print the command line used for running this program.
        logfile : str, optional
                  Path to log file.

        Returns
        -------
        outfile : str
                A VCF file.
        """

        if self.shapeit_folder is None:
            raise Exception("shapeit_folder must be defined")

        outfile = "{0}.vcf".format(output_prefix)

        Arg = namedtuple('Argument', 'option value')

        args = [Arg('--input-haps',
                    '{0}.gz {0}.sample'.format(input_prefix)), Arg('--output-vcf', outfile)]

        if logfile is not None:
            args.append(Arg('--output-log', logfile))

        compressRunner = None
        if compress is True:
            runner = RunProgram(path=self.shapeit_folder, program='shapeit -convert', args=args)
            if verbose is True: print("Command line is: {0}".format(runner.cmd_line))
            runner.run_checkoutput()
            compressRunner = RunProgram(path=self.bgzip_folder,
                                        program='bgzip',
                                        parameters=['-c', outfile, '>', outfile+".gz"])
            compressRunner.run_checkoutput()
            os.remove(outfile)
            outfile = outfile+".gz"
        else:
            runner = RunProgram(path=self.shapeit_folder,
                                program='shapeit -convert',
                                args=arguments)
            if verbose is True:
                print("Command line is: {0}".format(runner.cmd_line))
            runner.run_checkoutput()

        return outfile
