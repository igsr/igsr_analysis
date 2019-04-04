"""
Created on 21 March 2019

@author: galdam
"""

import eHive

import random
import string
import os
import shutil
import json
from typing import Dict, List, Tuple
from Utils.RunProgram import RunProgram


class Singularity(eHive.BaseRunnable):
    """
    The singularity runnable is intended as a base class for running workflow tasks inside singularity images.
    In most cases, it's sufficient to override the following parameters:
    - PIPELINE: The name of the pipeline module. Used in logging and as part of the default prefix
    - CMD: a new-style python format string. This will be the command that should be run in the singularity image.
    - CMD_ARGS: required arguments. These will be searched for in the 
    """
    PIPELINE = None
    CMD = None
    CMD_ARGS = []
    CMD_KWARGS = []
    FILES = dict()

    def run(self):
        """

        :return:
        """
        self.setup()

        # Setup the working directory
        self.open_working_dir()

        # Some logging
        if self.debug:
            with open(f"{self.working_dir}_log.txt", 'w') as log:
                log.write(str(self._BaseRunnable__params.__dict__))
                log.write('\n\n')
                log.write(json.dumps(self, default=lambda o: getattr(o, '__dict__', str(o))))
                log.write('\n\n')

        # Build the commands
        s_cmd = self.make_singularity_command()
        t_cmd = self.make_task_command()
        cmd = f"{s_cmd} {t_cmd}"

        # Some logging
        if self.debug:
            with open(f"{self.working_dir}_log.txt", 'a') as log:
                log.write(f"Running Command: {cmd}\n\n")

        # Execute the command
        self.execute_program(cmd)

        # Store the outputs and delete the temporary directory
        self.close_working_dir()

    # ---------------------------------------------------------------------------------------------------
    #
    def setup(self):
        """
        Setup the general parameters: define the output directory, prefix and verify the file list.
        """
        assert self.PIPELINE is not None, "Pipeline name has not been defined."

        self.output_directory = self.get_output_directory()
        self.working_dir = None
        self.prefix = self.get_prefix()
        self.get_output_file_list()

    def get_output_directory(self) -> str:
        """

        :return: str, the path to the output directory
        """
        root_output = self.param_required('root_output_dir')
        if self.param_exists('dir_label_params'):
            dir_label_params = self.param('dir_label_params')
            output_dir = os.path.join(root_output, *[self.param(d).replace(':', '.') for d in dir_label_params])
        else:
            output_dir = root_output
        return output_dir

    def get_prefix(self) -> str:
        """

        :return: str, the prefix for the
        """
        prefix = f"{self.param_required('basename')}.{self.PIPELINE}"
        prefix = prefix.replace(':', '.')
        return prefix

    def get_output_file_list(self) -> Dict[str, str]:
        """

        :return: {file label: file path, ...}
        """
        return {k: os.path.join(self.output_directory, v.format(PREFIX=self.prefix)) for k, v in self.FILES.items()}

    def make_singularity_command(self) -> str:
        """
        Sets up the singularity command. Includes selection of the executable, image and setting the working directory.
        :return: str, command for singularity
        """
        assert self.working_dir is not None, "Working directory has not been defined"
        command = [
            self.param_required('singularity_executable'),
            "exec",
            f"--pwd {self.working_dir}",
            f"{self.param_required('singularity_cache')}/{self.param_required('singularity_image')}",
        ]
        command = ' '.join(command)
        return command

    def make_task_command(self) -> str:
        """

        :return: str, command within singularity
        """
        cmd = self.CMD.format(WORKING_DIR=self.working_dir, PREFIX=self.prefix,
                              ARGS=self.unpack_cmd_kwargs(), **self.get_cmd_args())
        return cmd

    def get_cmd_args(self) -> Dict[str, str]:
        """

        :return:
        """
        return {k: self.param(k) for k in
                self.CMD_ARGS if self.param_is_defined(k)}

    def unpack_cmd_kwargs(self) -> str:
        """
        Utility for unpacking key-value parameters in the task command.
        :return: str, parameters
        """
        return ' '.join([f"--{k} {self.param(k)}" for k in self.CMD_KWARGS if self.param_is_defined(k)])

    def execute_program(self, cmd: str):
        """

        :param cmd:
        :return:
        """
        # Run the commands
        rp = RunProgram(cmd_line=cmd)
        rp.run_checkoutput()

    def open_working_dir(self):
        """
        Create and define the working directory
        :return:
        """
        working_dir = os.path.join(
            self.output_directory,
            f"{self.prefix}_{random_suffix()}"
        )
        os.makedirs(working_dir)
        self.working_dir = working_dir

    def close_working_dir(self) -> List[Tuple[str, str]]:
        """
        Move the files from the temporary directory to the output directory.
        This step attempts to be atomic; if the files can not be moved,
        all files will be returned to the temporary directory
        :return: list, moved files
        """
        files = os.listdir(self.working_dir)
        transactions = [(os.path.join(self.working_dir, file_name), os.path.join(self.output_directory, file_name))
                        for file_name in files]
        # Check first (it violates the pythonic principles of EAFP, but is cleaner)
        conflicts = [dst for src, dst in transactions if os.path.exists(dst)]
        if conflicts:
            shutil.Error("Cannot close working directory, files already exist: {}".format(
                ' ; '.join(conflicts)))

        for src, dst in transactions:
            try:
                shutil.move(src, dst)
            except shutil.Error as err:
                raise shutil.Error(f"Cannot close working directory, file already exist: {dst}") from err
        if not self.debug:
            shutil.rmtree(self.working_dir)
        return transactions

    def write_output(self):
        """
        Define thr output to pass back to eHive
        """
        self.warning('Work is done!')
        output_files = self.get_output_file_list()
        for f in output_files.values():
            assert os.path.exists(f), f"Missing file: {f}"
        self.dataflow(output_files, 1)


def random_suffix(suffix_length=5):
    """
    Generates a randomised string of numbers and lowercase characters.
    This can be used to differentiate temporary directories
    :param suffix_length: int, the number of characters to generate, default: 5
    :return: string
    """
    return ''.join(random.choices(string.ascii_lowercase + string.digits, k=suffix_length))

