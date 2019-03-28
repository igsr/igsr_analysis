"""
Created on 21 March 2019

@author: galdam
"""

import os
from Utils.RunSingularity import Singularity


class Test(Singularity):
    """
    Runs the Star alignment task in the broadinstitute/gtex_rnaseq singularity image.
    """
    CMD_ARGS = ['star_index', 'num_threads']
    CMD = "touch {WORKING_DIR}/$(basename {star_index}).touched"
    FILES = None
    PIPELINE = 'Test'

    def write_output(self):
        self.FILES = {'star_index': os.path.join(self.output_directory, 'INDEXAGAIN')}
        super().write_output()

    def execute_program(self, cmd):
        cmd = self.make_task_command()
        super().execute_program(cmd)






