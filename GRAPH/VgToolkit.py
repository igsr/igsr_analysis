'''
Created on 14 Sep 2021

@author: ernesto lowy ernestolowy@gmail.com
'''
import pdb
import glob

from Utils.RunProgram import RunProgram
from collections import namedtuple

class VG(object):
    """
    Class to run the different programs within the vg-toolkit
    (https://github.com/vgteam/vg)

    Class variables
    ---------------
    vg_folder : str, Optional
                Path to folder containing the vg binaries.
    """
    vg_folder = None 

    def __init__(self):
        pass
    
    def run_autoindex(self, ref: str, vcf: str, prefix: str, verbose: bool=False):
        """
        run vg autoindex

        Parameters
        ----------
        ref : str
              FASTA file containing the reference sequence
        vcf : str
              VCF file with sequence names matching -r
        prefix : str
                 Output prefix
        verbose : bool,  default=False
                  if true, then print the command line used for running this program

        Returns
        -------
        outfiles : list
                   List of output files
        """
        Arg = namedtuple('Argument', 'option value')
        args = [Arg('--workflow', 'giraffe')]
        args.append(Arg('-r', ref))
        args.append(Arg('-v', vcf))
        args.append(Arg('-p', prefix))

        program_cmd= f"{VG.vg_folder}/vg autoindex" if VG.vg_folder else "vg autoindex"
        vg_runner = RunProgram(program=program_cmd,
                               args=args)
        
        if verbose is True:
            print("Command line is: {0}".format(vg_runner.cmd_line))

        stdout, stderr, is_error = vg_runner.run_popen(raise_exc=False)

        outfiles = glob.glob(f"{prefix}*")
        
        return outfiles

    def run_giraffe(self, min: str, dist: str, fastq: str, prefix: str, gbz_f: str=None, gwbt_f: str=None, gbwt_g: str=None, verbose: bool= False ) -> str:
        """
        run vg giraffe

        Parameters
        ----------
        gbz_f : str, Optional
                path to GBZ file (GBWT index + GBWTGraph)
                Not required if gbwt_f and gbwt_g are provided
        gbwt_f : str, Optional.
                 path to GBWT index file. Not required if gbz_f is provided
        gbwt_g : str, Optional 
                 path to GBWTGraph file. Not required if gbz_f is provided
        min : str
             path to minimizer index file
        dist : str
              path to the distance index file
        fastq : str
                path to FASTQ file
        prefix : str
                 Output prefix
        verbose : bool,  default=False
                  if true, then print the command line used for running this program
        
        Returns
        -------
        gam file : str
                   Path to gam file
        """
        Arg = namedtuple('Argument', 'option value')
        args = (Arg('-Z', gbz_f), Arg('-m', min), Arg('-d', dist), Arg('-f', fastq), Arg('>', f"{prefix}.gam"))

        program_cmd= f"{VG.vg_folder}/vg giraffe" if VG.vg_folder else "vg giraffe"

        vg_runner = RunProgram(program=program_cmd,
                               args=args)
        
        if verbose is True:
            print("Command line is: {0}".format(vg_runner.cmd_line))

        stdout, stderr, is_error = vg_runner.run_popen(raise_exc=False)

        return f"{prefix}.gam"
 
    def run_stats(self, aln_f: str, verbose: bool= False ) -> str:
        """
        run vg stats

        Parameters
        ----------
        aln_f : str
               path to alignment file
        
        Returns
        -------
        stats_f : str
                  Path to file containing stats on the alingment
        """
        Arg = namedtuple('Argument', 'option value')
        args = (Arg('-a', aln_f), Arg('>', f"{aln_f}.stats"))

        program_cmd= f"{VG.vg_folder}/vg stats" if VG.vg_folder else "vg stats"

        vg_runner = RunProgram(program=program_cmd,
                               args=args)
        
        if verbose is True:
            print("Command line is: {0}".format(vg_runner.cmd_line))

        stdout, stderr, is_error = vg_runner.run_popen(raise_exc=False)

        return f"{aln_f}.stats"

    def run_augment(self, vg_f: str, aln_f: str, prefix: str, verbose: bool=False):
        """
        run vg augment

        Parameters
        ----------
        vg_f : str
               Path to graph.vg file
        aln_f : str
                Path to aln.gam file
        prefix : str
                 Output prefix
        
        Returns
        -------
        aug_graph_f : str
                      Path to augmented.vg file
        aug_aln_f : str
                    Path to augmented.gam file
        """
        Arg = namedtuple('Argument', 'option value')
        args = (Arg('', vg_f), Arg('', aln_f), Arg('-A', f"{prefix}.gam"), Arg('>', f"{prefix}.vg"))

        program_cmd= f"{VG.vg_folder}/vg augment" if VG.vg_folder else "vg augment"

        vg_runner = RunProgram(program=program_cmd,
                               args=args)
        
        if verbose is True:
            print("Command line is: {0}".format(vg_runner.cmd_line))

        stdout, stderr, is_error = vg_runner.run_popen(raise_exc=False)

        return f"{prefix}.vg", f"{prefix}.gam"
    
    def run_pack(self, vg_f: str, aln_f: str, prefix: str, verbose: bool=False, **kwargs) -> str:
        """
        run vg pack

        Parameters
        ----------
        vg_f : str
               Path to graph.vg file
        aln_f : str
                Path to aln.gam file
        prefix : str
                 Output prefix
        verbose : bool,  default=False
                  if true, then print the command line used for running this program
        **kwargs: Arbitrary keyword arguments.
        
        Returns
        -------
        {prefix}.pack : str
                      Path to .pack file
        """
        allowed_keys = ['Q'] # allowed arbitrary args

        Arg = namedtuple('Argument', 'option value')
        args = [Arg('-x', vg_f), Arg('-g', aln_f), Arg('-o', f"{prefix}.pack")]

        ## add now the **kwargs
        args.extend([Arg(f"-{k}", v) for k, v in kwargs.items() if k in allowed_keys])


        program_cmd= f"{VG.vg_folder}/vg pack" if VG.vg_folder else "vg pack"

        vg_runner = RunProgram(program=program_cmd,
                               args=args)
        
        if verbose is True:
            print("Command line is: {0}".format(vg_runner.cmd_line))

        stdout, stderr, is_error = vg_runner.run_popen(raise_exc=False)

        return f"{prefix}.pack"
    
    def run_call(self, vg_f: str, pack_f: str, prefix: str, verbose: bool=True) -> str:
        """
        run vg pack

        Parameters
        ----------
        vg_f : str
               Path to graph.vg file
        pack_f : str
                Path to aln.pack file
        prefix : str
                 Output prefix
        verbose : bool,  default=False
                  if true, then print the command line used for running this program
        
        Returns
        -------
        {prefix}.vcf : str
                      Path to .vcf file
        """
        allowed_keys = ['Q'] # allowed arbitrary args

        Arg = namedtuple('Argument', 'option value')
        args = [Arg('', vg_f), Arg('-k', pack_f), Arg('>', f"{prefix}.vcf")]

        program_cmd= f"{VG.vg_folder}/vg call" if VG.vg_folder else "vg call"

        vg_runner = RunProgram(program=program_cmd,
                               args=args)
        
        if verbose is True:
            print("Command line is: {0}".format(vg_runner.cmd_line))

        stdout, stderr, is_error = vg_runner.run_popen(raise_exc=False)

        return f"{prefix}.vcf"




