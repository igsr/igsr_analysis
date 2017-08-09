'''
Created on 03 Aug 2017

@author: ernesto
'''

import os
import subprocess
import pdb

class VcfUtils(object):
    '''
    Class to represent a misc of actions that can be done on a VCF file
    '''


    def __init__(self, vcf, bcftools_folder=None):
        '''
        Constructor

        Class variables
        ---------------
        vcf : str, Required
             Path to gzipped vcf file
        bcftools_folder : str, Required
                         Path to folder containing the bcftools binary
        '''

        if os.path.isfile(vcf) is False:
            raise Exception("File does not exist")

        self.vcf = vcf
        self.bcftools_folder = bcftools_folder

    def reheader(self, newheader, outprefix, samplefile=None, verbose=False):
        '''
        Modifiy the VCF's header with the newheader

        Parameters
        ----------
        newheader : string, required
                   Path to the file containing the new header
        outprefix : string, required
                    prefix for output files
        samplefile : string, optional
                     Path to the file with the sample names that will included
                     in the new header
        verbose : bool, optional
                  increase the verbosity, default=False

        Returns
        -------
        Path to the VCF with the modified header
        '''

        command = ""
        if self.bcftools_folder:
            command += self.bcftools_folder + "/"

        outfile=outprefix+".reheaded.vcf.gz"

        if samplefile is not None:
            command += "bcftools reheader -h {0} -o {1} -s {2} {3}".format(newheader, outfile, samplefile, self.vcf)
        else:
            command += "bcftools reheader -h {0} -o {1} {2}".format(newheader, outfile, self.vcf)

        if verbose is True:
            print("Command is: %s" % command)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            print("Cmd used was {0}".format(exc.cmd))
            raise Exception(exc.output)

        return outfile

