"""
Created on 21 March 2019

@author: galdam
"""
from typing import Dict
from Utils.RunSingularity import Singularity


class MarkDuplicates(Singularity):
    """
    Runs the Mark Duplicates task in the broadinstitute/gtex_rnaseq singularity image.
    """
    PIPELINE = 'markdups'
    CMD_ARGS = ['bam_file', 'memory']
    CMD = ("python3 -u /src/run_MarkDuplicates.py {bam_file} {PREFIX} -o {WORKING_DIR} --memory {memory}")
    FILES = {
        'md_bam_file': None,
        'metrics': "{PREFIX}.marked_dup_metrics.txt"
    }

    def get_output_file_list(self) -> Dict[str, str]:
        out_bam = self.param('bam_file')
        out_bam = out_bam.rsplit('.', 1)[0] + '.md.bam'
        self.FILES['md_bam_file'] = out_bam
        return super().get_output_file_list()





