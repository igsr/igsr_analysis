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
                Path to folder containing the vg binaries
    arg : namedtuple
          Containing a particular argument and its value
    """
    vg_folder = None 
    arg = namedtuple('Argument', 'option value')

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
        args = [VG.arg('--workflow', 'giraffe')]
        args.append(VG.arg('-r', ref))
        args.append(VG.arg('-v', vcf))
        args.append(VG.arg('-p', prefix))

        program_cmd= f"{VG.vg_folder}/vg autoindex" if VG.vg_folder else "vg autoindex"
        vg_runner = RunProgram(program=program_cmd,
                               args=args)
        
        if verbose is True:
            print("Command line is: {0}".format(vg_runner.cmd_line))

        stdout, stderr, is_error = vg_runner.run_popen(raise_exc=True)

        outfiles = glob.glob(f"{prefix}*")
        
        return outfiles

    def run_giraffe(self, fastq: str, prefix: str,  verbose: bool= False, **kwargs ) -> str:
        """
        run vg giraffe

        Parameters
        ----------
        fastq : str
                Path to FASTQ or comma-separated FASTQ files (if pair is provided) 
        prefix : str
                 Output prefix
        verbose : bool
                  if true, then print the command line used for running this program
        **kwargs: Arbitrary keyword arguments. Check the `vg giraffe` help for
                  more information.
        
        Returns
        -------
        gam file : str
                   Path to gam file
        """
        allowed_keys = ['H', 'Z', 'm', 'd', 'g', 't' ] # allowed arbitrary args

        files = fastq.split(",")
        files = [x.strip() for x in files]
        if len(files)>1:
            args = [VG.arg('-f', files[0]), VG.arg('-f', files[1])]
        else:
            args = [VG.arg('-f', files[0])]

        ## add now the **kwargs
        args.extend([VG.arg(f"-{k}", v) for k, v in kwargs.items() if k in allowed_keys])
        args.append(VG.arg('>', f"{prefix}.gam"))

        program_cmd= f"{VG.vg_folder}/vg giraffe" if VG.vg_folder else "vg giraffe"

        vg_runner = RunProgram(program=program_cmd,
                               args=args)
        
        if verbose is True:
            print("Command line is: {0}".format(vg_runner.cmd_line))

        stdout, stderr, is_error = vg_runner.run_popen(raise_exc=True)

        return f"{prefix}.gam"
 
    def run_chunk(self, *params, verbose: bool= True, **kwargs)->None:
        """
        Splits a graph and/or alignment into smaller chunks

        Parameters
        ----------
        verbose : bool
                  if true, then print the command line used for running this program.
        *params: Variable length parameter list.
        **kwargs: Arbitrary keyword arguments. Check the `vg chunk` help for
                  more information.
        
        Returns
        -------
        outfiles : list
                   Files having the different chunks
        """
        ## add the **kwargs
        allowed_keys = ['x', 'O', 'b', 'n'] # allowed arbitrary args
        args = ([VG.arg(f"-{k}", v) for k, v in kwargs.items() if k in allowed_keys])

        # mandatory parameters
        if not 'b' in kwargs:
            raise Exception("'b' argument is required by this method")        

        # add now the *params
        params_chunk = []
        allowed_opts = ['M'] # allowed options
        params_chunk.extend([f"-{x}" for x in params if x in allowed_opts])
        ext = f"{kwargs['O']}" if 'O' in kwargs else 'vg'

        program_cmd= f"{VG.vg_folder}/vg chunk" if VG.vg_folder else "vg chunk"
        vg_runner = RunProgram(program=program_cmd,
                               parameters=params_chunk,
                               args=args)
        
        if verbose is True:
            print("Command line is: {0}".format(vg_runner.cmd_line))

        stdout, stderr, is_error = vg_runner.run_popen(raise_exc=True)

        # get output files
        outfiles = glob.glob(f"{kwargs['b']}*.{ext}")
        return outfiles

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
        args = (VG.arg('-a', aln_f), VG.arg('>', f"{aln_f}.stats"))

        program_cmd= f"{VG.vg_folder}/vg stats" if VG.vg_folder else "vg stats"

        vg_runner = RunProgram(program=program_cmd,
                               args=args)
        
        if verbose is True:
            print("Command line is: {0}".format(vg_runner.cmd_line))

        stdout, stderr, is_error = vg_runner.run_popen(raise_exc=True)

        return f"{aln_f}.stats"

    def run_augment(self, *params, graph_f: str, aln_f: str, prefix: str, verbose: bool=True, 
                    **kwargs):
        """
        run vg augment

        Parameters
        ----------
        graph_f : str
                  Path to graph file
        aln_f : str
                Path to aln.gam file
        prefix : str
                 Output prefix
        *params: Variable length parameter list.
        **kwargs: Arbitrary keyword arguments. Check the `vg augment` help for
                  more information.
        
        Returns
        -------
        aug_graph_f : str
                      Path to {prefix}.aug.pg file
        aug_aln_f : str
                    Path to {prefix}.aug.gam
        """
        ## add the **kwargs
        allowed_keys = ['Q', 'q', 'm'] # allowed arbitrary args
        args = ([VG.arg(f"-{k}", v) for k, v in kwargs.items() if k in allowed_keys])

        # add the mandatory argument
        args.append(VG.arg('-A', f"{prefix}.aug.gam"))

        # add the mandatory params
        params_aug = [graph_f, aln_f]

        # add now the *params
        allowed_opts = ['s'] # allowed options
        params_aug.extend([f"-{x}" for x in params if x in allowed_opts])
        params_aug.extend(['>', f"{prefix}.aug.pg"])


        program_cmd= f"{VG.vg_folder}/vg augment" if VG.vg_folder else "vg augment"

        vg_runner = RunProgram(program=program_cmd,
                               args=args,
                               parameters=params_aug)
        
        if verbose is True:
            print("Command line is: {0}".format(vg_runner.cmd_line))

        stdout, stderr, is_error = vg_runner.run_popen(raise_exc=True)

        return f"{prefix}.aug.pg", f"{prefix}.aug.gam"
    
    def run_pack(self, vg_f: str, aln_f: str, prefix: str, verbose: bool=True, **kwargs) -> str:
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
        **kwargs: Arbitrary keyword arguments. Check the `vg pack` help for
                  for more information.
        
        Returns
        -------
        {prefix}.pack : str
                      Path to .pack file
        """
        allowed_keys = ['Q'] # allowed arbitrary args

        args = [VG.arg('-x', vg_f), VG.arg('-g', aln_f), VG.arg('-o', f"{prefix}.pack")]

        ## add now the **kwargs
        args.extend([VG.arg(f"-{k}", v) for k, v in kwargs.items() if k in allowed_keys])

        program_cmd= f"{VG.vg_folder}/vg pack" if VG.vg_folder else "vg pack"

        vg_runner = RunProgram(program=program_cmd,
                               args=args)
        
        if verbose is True:
            print("Command line is: {0}".format(vg_runner.cmd_line))

        stdout, stderr, is_error = vg_runner.run_popen(raise_exc=True)

        return f"{prefix}.pack"
    
    def run_call(self, graph_f: str, prefix: str, verbose: bool=False, **kwargs) -> str:
        """
        run vg pack

        Parameters
        ----------
        graph_f : str
                  Path to graph file
        prefix : str
                 Output prefix
        verbose : bool,  default=False
                  if true, then print the command line used for running this program.
        **kwargs: Arbitrary keyword arguments. Check the `vg pack` help for
                  for more information.
        
        Returns
        -------
        {prefix}.vcf : str
                      Path to .vcf file
        """
        allowed_keys = ['k', 'r', 's'] # allowed arbitrary args

        ## add now the **kwargs
        args = [VG.arg(f"-{k}", v) for k, v in kwargs.items() if k in allowed_keys]

        parameters = [graph_f, '>', f"{prefix}.vcf"]
        
        program_cmd= f"{VG.vg_folder}/vg call" if VG.vg_folder else "vg call"

        vg_runner = RunProgram(program=program_cmd,
                               args=args,
                               parameters = parameters)
        
        if verbose is True:
            print("Command line is: {0}".format(vg_runner.cmd_line))

        stdout, stderr, is_error = vg_runner.run_popen(raise_exc=True)

        return f"{prefix}.vcf"

    def run_snarls(self, graph_f: str, prefix: str, verbose: bool=True):
        """
        run vg snarls

        Parameters
        ----------
        graph_f : str
                  Path to Graph file
        prefix : str
                 Output prefix
        verbose : bool,  default=False
                  if true, then print the command line used for running this program
        
        Returns
        -------
        {prefix}.snarls : str
                          Path to .snarls file
        """
        parameters = [graph_f, '>', f"{prefix}.snarls"]

        program_cmd= f"{VG.vg_folder}/vg snarls" if VG.vg_folder else "vg snarls"

        vg_runner = RunProgram(program=program_cmd,
                               parameters=parameters)
        
        if verbose is True:
            print("Command line is: {0}".format(vg_runner.cmd_line))

        stdout, stderr, is_error = vg_runner.run_popen(raise_exc=True)

        return f"{prefix}.snarls"