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

    def run_shapeit(self, input_gen, input_init, input_scaffold, input_map, output_prefix, from=None, to=None, **kwargs):
        '''
        Run Shapeit
        
        Parameters
        ----------
        input_gen : Tuple, Required
                    specifies the genotype/GL input data that you obtain from Beagle4, i.e. ('input.shapeit.20.gen.gz','input.shapeit.20.gen.sample')
        input_init : Tuple, Required
                     specifies the haplotypes that you obtain from Beagle4, i.e. ('input.shapeit.20.hap.gz','input.shapeit.20.hap.sample')
        input_scaffold : Tuple, Required
                         SNP-array derived haplotype scaffold used by SHAPEIT. It has to be in Impute2 format. i.e. ('scaffold.haps.gz','scaffold.haps.sample')
        input_map : str, Required
                    Path with the genetic map
        output_prefix : str, Required
                        Prefix used for the 2 output files estimated by SHAPEIT, i.e. 'output.shapeit.20.haps.gz output.shapeit.20.haps.sample',
        from : int, Optional
               specify the region to be phased
        to : int, Optional
             specify the region to be phased
        
        '''

        command = ""

        if self.shapeit_folder:
            command += self.shapeit_folder+"/"
