'''
Created on 21 Jul 2017

@author: ernesto
'''
import os
import pdb
import subprocess

class Shapeit(object):
    '''
    Class to run SHAPEIT
    '''

    def __init__(self, shapeit_folder=None):
        '''
        Constructor

        Class variables
        ---------------
        shapeit_folder : str, Optional
                         Path to folder containing Shapeit binary
        '''
        self.shapeit_folder = shapeit_folder

    def run_shapeit(self, input_gen, input_init, input_scaffold, output_prefix, input_map=None,
                    verbose=False, **kwargs):
        '''
        Run Shapeit
        
        Parameters
        ----------
        input_gen : str, Required
                    specifies the genotype/GL input data that you obtain from Beagle4, i.e. 'input.shapeit.20.gen.gz input.shapeit.20.gen.sample'
        input_init : str, Required
                     specifies the haplotypes that you obtain from Beagle4, i.e. 'input.shapeit.20.hap.gz input.shapeit.20.hap.sample'
        input_scaffold : str, Required
                         SNP-array derived haplotype scaffold used by SHAPEIT. It has to be in Impute2 format. i.e. 'scaffold.haps.gz scaffold.haps.sample'
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
        
        '''

        

        command = ""

        if self.shapeit_folder:
            command += self.shapeit_folder+"/"

        command += "shapeit -call --input-gen {0} --input-init {1} --input-scaffold {2} --output-max {3}.haps.gz {3}.haps.sample --output-log {3}.log".format(input_gen, input_init, input_scaffold,
                                                                                                                                                              output_prefix)
        for k,v in kwargs.items():
            command += " --{0} {1}".format(k,v)

        if verbose==True:
            print("Command used was: %s" % command)
             
        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            print("Command used was: {0}".format(command))
            raise Exception(exc.output)
