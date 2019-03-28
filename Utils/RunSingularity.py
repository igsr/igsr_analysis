"""
Created on 21 March 2019

@author: galdam
"""

import eHive

import random
import string
import os
import shutil
from Utils.RunProgram import RunProgram
import json


class Singularity(eHive.BaseRunnable):
    CMD_KWARGS = []
    CMD_ARGS = []
    CMD = None
    FILES = dict()
    PIPELINE = None

    def run(self):
        self.debug = True
        self.output_directory = self.get_output_directory()
        self.working_dir = None
        self.prefix = self.get_prefix()

        # Setup the working directory
        self.open_working_dir()

        # Some logging
        if self.debug:
            with open(f"{self.working_dir}_log.txt", 'w') as log:
                try:
                    log.write(str(self._BaseRunnable__params.__dict__))
                except:
                    pass
                log.write('\n\n')
                log.write(json.dumps(self, default=lambda o: getattr(o, '__dict__', str(o))))
                log.write('\n\n')


        # Build the commands
        s_cmd = self.make_singularity_command()
        t_cmd = self.make_task_command()
        cmd = f"{s_cmd} {t_cmd}"

        # Some logging
        if self.debug:
            with open(f"{self.working_dir}_log.txt", 'w') as log:
                try:
                    log.write(str(self._BaseRunnable__params.__dict__))
                except:
                    pass
                log.write('\n\n')
                log.write(json.dumps(self, default=lambda o: getattr(o, '__dict__', str(o))))
                log.write('\n\n')
                log.write(f"Running Command: {cmd}\n\n")

        self.execute_program(cmd)

        # Store the outputs and delete the temporary directory
        self.close_working_dir()

    def execute_program(self, cmd):
        # Run the commands
        rp = RunProgram(cmd_line=cmd)
        rp.run_checkoutput()

    def make_singularity_command(self):
        """
        Sets up the singularity command. Includes selection of the executable, image and setting the working directory.
        :return: string, command for singularity
        """
        assert self.working_dir is not None, "Working directory has not been defined"
        command = [
            self.param_required('singularity_executable'),
            "exec",
            f"{self.param_required('singularity_cache')}/{self.param_required('singularity_image')}",
            f"--pwd {self.working_dir}"
        ]
        command = ' '.join(command)
        return command

    def make_task_command(self):
        """
        Override this to define a task to run through singularity
        """
        cmd = self.CMD.format(WORKING_DIR=self.working_dir, PREFIX=self.prefix,
                              ARGS=self.unpack_cmd_kwargs(), **self.get_cmd_args())
        return cmd

    def get_cmd_args(self):
        return {k: self.param(k) for k in
                self.CMD_ARGS if self.param_is_defined(k)}

    def unpack_cmd_kwargs(self):
        """
        Utility for unpacking key-value prarmeters in the task command.
        :return: String, parameters
        """
        return ' '.join([f"--{k} {self.param(k)}" for k in self.CMD_KWARGS if self.param_is_defined(k)])

    # ---------------------------------------------------------------------------------------------------
    #
    def get_prefix(self):
        return f"{self.param_required('basename')}.{self.PIPELINE}"

    def get_output_directory(self):
        root_output = self.param_required('root_output_dir')
        if self.param_exists('dir_label_params'):
            dir_label_params = self.param('dir_label_params')
            output_dir = os.path.join(root_output, *[self.param(d) for d in dir_label_params])
        else:
            output_dir = root_output
        return output_dir

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

    def close_working_dir(self):
        """
        Move the files from the temporary directory to the output directory.
        This step attempts to be atomic; if the files can not be moved,
        all files will be returned to the temporary directory
        :return: list, moved files
        """
        files = os.listdir(self.working_dir)
        target_files = os.listdir(self.output_directory)
        if set(files) & set(target_files):
            shutil.Error("Cannot close working directory, files already exist: {}".format(
                ' ; '.join(set(files) & set(target_files))))
        moved_files = []

        try:
            for file_name in files:
                shutil.move(os.path.join(self.working_dir, file_name), self.output_directory)
                moved_files.append(os.path.join(self.output_directory, file_name))
        # If one of the files already exists, move the files back before raising the original error.
        # Use print errors for any new errors encountered at this point
        except shutil.Error as err:
            for file_path in moved_files:
                try:
                    shutil.move(file_path, self.output_directory)
                except shutil.Error as other_err:
                    print(f"Error raised whilst solving file conflict: {other_err}")
            raise shutil.Error("Cannot close working directory, files already exist") from err
        if not self.debug:
            shutil.rmtree(self.working_dir)
        return moved_files

    def write_output(self):
        self.warning('Work is done!')
        files = {k: os.path.join(self.output_directory, v.format(self.prefix)) for k, v in self.FILES.items()}
        self.dataflow(files, 1)


def random_suffix(suffix_length=5):
    """
    Generates a randomised string of numbers and lowercase characters.
    This can be used to differentiate temporary directories
    :param suffix_length: int, the number of characters to generate, default: 5
    :return: string
    """
    return ''.join(random.choices(string.ascii_lowercase + string.digits, k=suffix_length))

