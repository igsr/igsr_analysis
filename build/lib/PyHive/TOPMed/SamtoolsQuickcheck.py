"""
Created on 21 March 2019

@author: galdam
"""

from Utils.RunSingularity import Singularity


class SamtoolsQuickcheck(Singularity):
    """
    Runs the samtools quickcheck task in the broadinstitute/gtex_rnaseq singularity image.
    """
    PIPELINE = 'samtools_quickcheck'
    CMD_ARGS = ['bam_file_to_check']
    CMD = "samtools quickcheck {bam_file_to_check}"

