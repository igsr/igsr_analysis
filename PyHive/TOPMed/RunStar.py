"""
Created on 21 March 2019

@author: galdam
"""

import eHive
from Utils.RunSingularity import Singularity
#from Utils.TopMedRunnable import TopMedRunnable


class


class RunStar(eHive.BaseRunnable):


    def write_output(self):
        self.warning('Work is done!')
        outdict = self.param('outdict')



class Star(Singularity):
    """
    Runs the Star alignment task in the broadinstitute/gtex_rnaseq singularity image.
    """

    CMD_KWARGS = [
        "outFilterMultimapNmax", "alignSJoverhangMin", "alignSJDBoverhangMin", "outFilterMismatchNmax",
        "outFilterMismatchNoverLmax", "alignIntronMin", "alignIntronMax", "alignMatesGapMax", "outFilterType",
        "outFilterScoreMinOverLread", "outFilterMatchNminOverLread", "limitSjdbInsertNsj", "outSAMstrandField",
        "outFilterIntronMotifs", "alignSoftClipAtReferenceEnds", "quantMode", "outSAMattrRGline", "outSAMattributes",
        "chimSegmentMin", "chimJunctionOverhangMin", "chimOutType", "chimMainSegmentMultNmax", "sjdbFileChrStartEnd",
    ]

    def make_task_command(self):
        cmd = (
            "/src/run_STAR.py {star_index} {fastq1} {fastq2} {PREFIX} "
            "--output_dir {WORKING_DIR} --threads {num_threads} {ARGS}"
        ).format(WORKING_DIR=self.working_dir, PREFIX=self.prefix,
                 ARGS=self.unpack_cmd_arguments(), **self.cmd_arguments)
        return cmd
