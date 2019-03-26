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


DEFAULT = {
    'singularity_executable': 'singularity'
}


class Singularity(eHive.BaseRunnable):
    CMD_KWARGS = []

    def run(self):

        options_dict = {k: self.param(k) for k in self.CMD_KWARGS if self.param_is_defined(k)}

        fastq1, fastq2 = sorted(self.param('fastq'))
        options_dict['fastq1'] = fastq1
        options_dict['fastq2'] = fastq2
        options_dict['star_index'] = self.param_required('star_index')
        options_dict['num_threads'] = self.param_required('num_threads')

         self.setup(
            s_cache=self.param_required('singularity_cache'),
            s_image=self.param_required('singularity_image'),
            s_executable=self.param_required('singularity_exe'),
            output_directory=self.param_required('output_directory'),
            prefix=self.param_required('prefix'),
            cmd_arguments=options_dict)
        self.run()


    def setup(self, s_cache, s_image, output_directory, prefix, s_executable=None, cmd_arguments=None):
        """

        :param s_cache: the singularity cache location
        :param s_image: the singularity image name
        :param s_executable: the singularity executable (default: 'singularity'
        """
        self.singularity_cache = s_cache
        self.singularity_image = s_image
        self.singularity_executable = DEFAULT['singularity_executable'] \
            if s_executable is None else s_executable
        self.cmd_arguments = dict() if cmd_arguments is None else cmd_arguments

        self.output_directory = output_directory
        self.working_dir = None
        self.prefix = prefix
        self.debug = True

    def execute(self):
        """

        :return:
        """
        # Setup the working directory
        self.open_working_dir()

        # Build the commands
        s_cmd = self.make_singularity_command()
        t_cmd = self.make_task_command()
        cmd = f"{s_cmd} {t_cmd}"
        print(f"Running Command: {cmd}")

        # Run the commands
        rp = RunProgram(cmd_line=cmd)
        rp.run_checkoutput()

        # Store the outputs and delete the temporary directory
        self.close_working_dir()

    def make_singularity_command(self):
        """
        Sets up the singularity command. Includes selection of the executable, image and setting the working directory.
        :return: string, command for singularity
        """
        assert self.working_dir is not None, "Working directory has not been defined"
        command = [
            self.singularity_executable,
            "exec",
            f"{self.singularity_cache}/{self.singularity_image}",
            f"--pwd {self.working_dir}"
        ]
        command = ' '.join(command)
        return command

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

    def make_task_command(self):
        """
        Override this to define a task to run through singularity
        """
        raise NotImplementedError()

    def unpack_cmd_arguments(self):
        """
        Utility for unpacking key-value prarmeters in the task command.
        :return: String, parameters
        """
        return ' '.join([f"--{k} {v}" for k, v in self.cmd_arguments.items() if k in self.CMD_KWARGS])


def random_suffix(suffix_length=5):
    """
    Generates a randomised string of numbers and lowercase characters.
    This can be used to differentiate temporary directories
    :param suffix_length: int, the number of characters to generate, default: 5
    :return: string
    """
    return ''.join(random.choices(string.ascii_lowercase + string.digits, k=suffix_length))

