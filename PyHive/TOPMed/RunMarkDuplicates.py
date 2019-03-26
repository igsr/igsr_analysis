"""
Created on 21 March 2019

@author: galdam
"""

from Utils.RunSingularity import Singularity


class MarkDuplicates(Singularity):
    """
    Runs the Mark Duplicates task in the broadinstitute/gtex_rnaseq singularity image.
    """

    CMD_KWARGS = []

    def make_task_command(self):
        """

        :return:
        """
        cmd = ("python3 -u /src/run_MarkDuplicates.py "
               "{input_bam} {PREFIX} --memory {memory}"
               ).format(PREFIX=self.prefix, **self.cmd_arguments)
        return cmd



