'''
Created on 21 Jul 2017

@author: ernesto
'''
import os
import pdb
import subprocess
from Utils.RunProgram import RunProgram
from collections import namedtuple

class Shapeit(object):
    '''
    Class to run SHAPEIT
    '''

    def __init__(self, bgzip_folder=None, shapeit_folder=None, ligateHAPLOTYPES_folder=None):
        '''
        Constructor

        Class variables
        ---------------
        bgzip_folder : str, Optional
                       Path to folder containing the bgzip binary
        shapeit_folder : str, Optional
                         Path to folder containing the Shapeit binary
        ligateHAPLOTYPES_folder : str, Optional
                                  Path to folder containing the ligateHAPLOTYPES binary
        '''
        self.bgzip_folder = bgzip_folder
        self.shapeit_folder = shapeit_folder
        self.ligateHAPLOTYPES_folder = ligateHAPLOTYPES_folder

    def run_shapeit(self, output_prefix, input_gen=None, input_init=None, input_scaffold=None, input_bed=None, duohmm=False, input_map=None,
                    verbose=False, **kwargs):
        '''
        Run Shapeit
        
        Parameters
        ----------
        input_gen : str, Optional
                    specifies the genotype/GL input data that you obtain from Beagle4, i.e. 'input.shapeit.20.gen.gz input.shapeit.20.gen.sample'
        input_init : str, Optional
                     specifies the haplotypes that you obtain from Beagle4, i.e. 'input.shapeit.20.hap.gz input.shapeit.20.hap.sample'
        input_scaffold : str, Optional
                         SNP-array derived haplotype scaffold used by SHAPEIT. It has to be in Impute2 format. i.e. 'scaffold.haps.gz scaffold.haps.sample'
        input_bed : str, Optional
                    Unphased genotypes in Plink Binary BED/BIM/FAM format. i.e. 'file.bed file.bim file.fam'
        duohmm : bool, Optional
                 If true, then activate the --duohmm option. Default: False
        output_prefix : str, Required
                        Prefix used for the 2 output files estimated by SHAPEIT, i.e. 'output.shapeit.20.haps.gz output.shapeit.20.haps.sample'
        input_map : str, Optional
                    Path with the genetic map
        i_from : int, Optional
               specify the region to be phased
        i_to : int, Optional
             specify the region to be phased
        verbose : bool, optional
                  if true, then print the command line used for running this program
        
        Returns
        -------
        A dict with the path to the 2 output files (*.haps.gz and *.haps.sample) that can be used with SHAPEIT
        '''

        if input_gen is None and input_bed is None:
            raise Exception("Error! Either --input-gen or --input-bed need to be specified as input for SHAPEIT")
            
        Arg = namedtuple('Argument', 'option value')

        args=[]

        if input_gen is not None:
            args.append(Arg('-call --input-gen',input_gen))
        elif input_bed is not None:
            args.append(Arg('--input-bed',input_bed))

        if input_init is not None:
            args.append(Arg('--input-init',input_init))

        if input_scaffold is not None:
            args.append(Arg('--input-scaffold',input_scaffold))

        if input_map is not None:
            args.append(Arg('--input-map',input_map))
        
        for k,v in kwargs.items():
            args.append(Arg('--{0}'.format(k),v))

        args.extend([Arg('--output-max','{0}.haps.gz {0}.haps.sample'.format(output_prefix)),Arg('--output-log','{0}.log'.format(output_prefix))])

        params=[]
        if duohmm is True: params=['--duohmm']

        runner=RunProgram(path=self.shapeit_folder, program='shapeit', args=args, parameters=params)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout=runner.run_checkoutput()
             
        outdict={ 
            'hap_gz' : '{0}.haps.gz'.format(output_prefix),
            'hap_sample' : '{0}.haps.sample'.format(output_prefix)
        }

        return outdict

    def ligate_shapeitchunks(self,vcf_f,scaffolded_samples,chunk_str,output_prefix,verbose=False):
        '''
        Run ligateHAPLOTYPES to ligate together all haplotype chunks produced by SHAPEIT
        (see https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#haplegsample)

        Parameters
        ----------
        vcf_f : str, Required
                VCF with the Genotype likelihoods
        scaffolded_samples : str, Required
                             File with the list of samples (separated by '\n') that have been scaffolded
        chunk_str : str, Required
                    String with the paths to the different files generated by SHAPEIT for the different 
                    chromosome chunks (i.e. 's2.chunk1.hap.gz s2.chunk1.hap.gz s2.chunk1.hap.gz')
        output_prefix : str, Required
                        String with the output prefixes (i.e. 'output.shapeit.22.ligated.haps.gz output.shapeit.22.ligated.haps.sample')

        Returns
        -------
        A dict with the path to the 2 output files (*.haps.gz and *.haps.sample)
        '''
        
        command = ""

        if self.ligateHAPLOTYPES_folder:
            command += self.ligateHAPLOTYPES_folder+"/"

        command += "ligateHAPLOTYPES --vcf {0} --scaffold {1} --chunks {2} --output {3}.ligated.haps.gz {3}.ligated.haps.sample".format(vcf_f,
                                                                                                                                        scaffolded_samples, 
                                                                                                                                        chunk_str,
                                                                                                                                        output_prefix)
        if verbose==True:
            print("Command used was: %s" % command)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            print("Command used was: {0}".format(command))
            raise Exception(exc.output)

        outdict={
            'hap_gz' : '{0}.ligated.haps.gz'.format(output_prefix),
            'hap_sample' : '{0}.ligated.haps.sample'.format(output_prefix)
        }

        return outdict

    def convert2vcf(self, input_prefix, output_prefix, compress=False, verbose=False, logfile=None):
        '''
        Function to use SHAPEIT's -convert in order to convert the *.haps.gz & *.haps.sample files into VCF

        Parameters
        ----------
        input_prefix : str, Required
                       Prefix for the files in HAPS/SAMPLE format
        output_prefix : str, Required
                        String with the output prefix for the VCF file
        verbose : bool, optional
                  if true, then print the command line used for running this program
        logfile : str, optional
                  Path for log file
        
        Returns
        -------
        A VCF file
        '''

        command = ""

        if self.shapeit_folder:
            command += self.shapeit_folder+"/"

        outfile="{0}.vcf".format(output_prefix)

        command += "shapeit -convert --input-haps {0}.gz {0}.sample --output-vcf {1}".format(input_prefix, outfile)

        if logfile is not None:
            command += " --output-log {0}".format(logfile)

        if verbose==True:
            print("Command used was: %s" % command)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            print("Command used was: {0}".format(command))
            raise Exception(exc.output)

        if compress is True:
            bgzip_path = ""
            if self.bgzip_folder:
                bgzip_path = "{0}/bgzip".format(self.bgzip_folder)
            else:
                bgzip_path = "bgzip"
            compress_cmd="{0} -c {1} > {1}.gz".format(bgzip_path, outfile)
            
            try:
                subprocess.check_output(compress_cmd, shell=True)
            except subprocess.CalledProcessError as exc:
                print("Command used was: {0}".format(compress_cmd))
                raise Exception(exc.output)
            
            os.remove(outfile)
            outfile = outfile+".gz"

        return outfile
