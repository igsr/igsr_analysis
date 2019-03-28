"""
Created on 21 March 2019

@author: galdam
"""

from Utils.RunSingularity import Singularity


class MarkDuplicates(Singularity):
    """
    Runs the Mark Duplicates task in the broadinstitute/gtex_rnaseq singularity image.
    """

    CMD_ARGS = ['input_bam', 'memory']
    CMD = ("python3 -u /src/run_MarkDuplicates.py {input_bam} {PREFIX} --memory {memory}")




