"""
Created on 21 March 2019

@author: galdam
"""

from Utils.RunSingularity import Singularity


class Rsem(Singularity):
    """
    Runs the Rsem task in the broadinstitute/gtex_rnaseq singularity image.
    """

    CMD_KWARGS = ['max_frag_len', 'estimate_rspd', 'is_stranded', 'paired_end']
    CMD_ARGS = ['num_threads', 'rsem_reference', 'transcriptome_bam']
    CMD = ("/src/run_RSEM.py {ARGS} --threads {num_threads} "
           "{rsem_reference} {transcriptome_bam} {PREFIX}")


