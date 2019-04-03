"""
Created on 21 March 2019

@author: galdam
"""

from Utils.RunSingularity import Singularity


class IndexBam(Singularity):
    """
    Runs the Rsem task in the broadinstitute/gtex_rnaseq singularity image.
    """
    PIPELINE = 'indexbam'
    CMD_ARGS = ['unindexed_bam_file']
    CMD = "samtools index {unindexed_bam_file}"

