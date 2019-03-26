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

    def make_task_command(self):
        """

        :return:
        """
        cmd = (
            "/src/run_RSEM.py {ARGS} --threads {num_threads} "
            "{rsem_reference} {transcriptome_bam} {PREFIX}"
        ).format(PREFIX=self.prefix, ARGS=self.unpack_cmd_arguments(),
                 **self.cmd_arguments)
        return cmd

