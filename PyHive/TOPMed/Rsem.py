"""
Created on 21 March 2019

@author: galdam
"""

from Utils.RunSingularity import Singularity


class Rsem(Singularity):
    """
    Runs the Rsem task in the broadinstitute/gtex_rnaseq singularity image.
    """
    PIPELINE = 'rsem'
    CMD_KWARGS = ['max_frag_len', 'estimate_rspd', 'is_stranded', 'paired_end']
    CMD_ARGS = ['num_threads', 'rsem_reference', 'transcriptome_bam']
    CMD = ("/src/run_RSEM.py {ARGS} -o {WORKING_DIR} --threads {num_threads} "
           "{rsem_reference} {transcriptome_bam} {PREFIX}")
    FILES = {
        'genes': "{PREFIX}.rsem.genes.results",
        'isoforms': "{PREFIX}.rsem.isoforms.results"
    }

    @property
    def prefix(self) -> str:
        """
        Overriding the function as the default output is in the format "<prefix>.rsem.*" and
        would otherwise result in <prefix>.rsem.rsem.*
        :return: str, the prefix for output files
        """
        return self.base_prefix
