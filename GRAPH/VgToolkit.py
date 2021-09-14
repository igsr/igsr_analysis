'''
Created on 14 Sep 2021

@author: ernesto lowy ernestolowy@gmail.com
'''
import pdb

from Utils.RunProgram import RunProgram

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

    def __init__(self, bam):
        """
        Constructor

        Parameters
        ----------
        bam : str
             Path to BAM file used for the variant calling process.
        """

        if os.path.isfile(bam) is False:
            raise Exception("BAM file does not exist")

        self.bam = bam
     
    def run_giraffe(self):

        program_cmd= f"{VG.vg_folder}/vg giraffe" if VG.vg_folder else "vg giraffe"

        runner = RunProgram(program=program_cmd,
                            parameters=['-h'])
