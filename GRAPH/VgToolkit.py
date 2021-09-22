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

    def run_giraffe(self, gbz_f: str, min: str, dist: str, fastq: str, prefix: str, verbose: bool= False ) -> str:
        """
        run vg giraffe

        Parameters
        ----------
        gbz_f : str
               path to GBZ file (GBWT index + GBWTGraph)
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

        pdb.set_trace()
        vg_runner = RunProgram(program=program_cmd,
                               args=args)
        
        if verbose is True:
            print("Command line is: {0}".format(vg_runner.cmd_line))

        stdout, stderr, is_error = vg_runner.run_popen(raise_exc=False)

        return f"{aln_f}.stats"

    def run_augment(self, vg_f: str, aln_f: str, prefix: str, verbose: bool=True):
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



