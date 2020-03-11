'''
Created on 03 May 2018

@author: ernesto
'''

import subprocess
import re
import os
from configparser import ConfigParser


class RunProgram(object):
    """
    SuperClass used to run an external program within a Python script

    Class Variables
    ---------------
    settings : ConfigParser, Required
               Object with configuration settings
               relevant for this class
    """

    def __init__(self, program, settingf, use_docker=False, path=None, args=None,
                 arg_sep=None, parameters=None, cmd_line=None, downpipe=None):
        """
        Constructor

        Instance Variables
        ------------------
        program : str, Required
                  Program to be run
        settingf : str, Required
                   Path to .ini file with settings
        use_docker : Bool, Optional
                     Use docker containter for running this command. Default: False
        path : str, optional
               Folder containing the 'program'
        args : list, optional
               List of named tuples tuples formed by an argument (or option) and
               its respective value: i.e. arg=namedtuple('Argument', 'option value')
        parameters : list, optional
                     List of parameters ['c','d']
        arg_sep : char, optional
                  char used as a separator between argument and value.
                  i.e. if '=' then we will get 'a'=1
                  Default is a single whitespace (i.e. ' ')
        cmd_line : str, optional
                   String with command line to run
        downpipe: list, optional
                  List of RunProgram objects that will be executed in a pipe after
                  self.program has been executed
        """
        pdb.set_trace()
        self.cmd_line = cmd_line
        self.program = program
        self.use_docker = use_docker
        if self.program is None and self.cmd_line is None:
            raise ValueError("Parameter 'cmd_line' or 'program' must be provided.")

        self.path = path
        self.args = args
        self.arg_sep = arg_sep if arg_sep is not None else ' '
        self.parameters = parameters
        self.downpipe = downpipe

        # parse settings file (in .ini file)
        parser = ConfigParser(allow_no_value=True)
        parser.optionxform = str

        parser.read(settingf)
        self.settings = parser

        # create the command line if is None
        if self.cmd_line is None:
            self.cmd_line = self.create_command_line()

    def create_command_line(self):
        """
        :returns
        str, cmd line
        """
        assert self.path is not None and self.use_docker is True, "Conflicting instructions, " \
                                                                  "do not know if use local" \
                                                                  " dependency or container"
        if self.path is not None:
            cmd_line = [os.path.join(self.path, self.program)]
        else:
            cmd_line = [self.program]

        # construct the command line
        if self.args is not None:
            cmd_line.append(' '.join([f"{option}{self.arg_sep}{parameter}"
                                      for option, parameter in self.args]))

        if self.parameters is not None:
            cmd_line.append(' '.join(["{0}".format(param) for param in self.parameters]))

        if self.downpipe is not None:
            cmd_line.append(' '.join(["| {0}".format(runO.cmd_line) for runO in self.downpipe]))

        return ' '.join(cmd_line)


    def run_popen(self, raise_exc=True):
        """
        Run self.program using subprocess Popen method
        (see https://docs.python.org/2/library/subprocess.html#module-subprocess)

        Parameters
        ----------
        raise_exc: bool, optional
                   If true, then raise an Exception when error is found. Default= True

        Returns
        -------
        tuple
             A tuple containing the STDOUT and STDERR from this program
        """
        log_f = None
        if self.settings.has_option('run_program', 'log_file'):
            log_f = open(self.settings.get('run_program', 'log_file'), 'w')

        # execute cmd_line
        p = subprocess.Popen(self.cmd_line, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, bufsize=256,
                             universal_newlines=True)

        # stderr
        patt = re.compile('#* ERROR|Error')
        is_exception = False
        stderr = ""
        for line in p.stderr:
            line = str(line.rstrip())
            stderr += line + "\n"
            if log_f is not None:
                log_f.write(line + "\n")
            m = patt.match(line)
            if m:
                is_exception = True

        # stdout
        stdout = ""
        for line in p.stdout:
            line = str(line.rstrip())
            stdout += line + "\n"

        if is_exception is True and raise_exc is True:
            raise Exception(stderr)

        if log_f is not None:
            log_f.close()

        return (stdout, stderr, is_exception)

    def run_checkoutput(self):
        """
        Run self.program using subprocess check_output method
        (see https://docs.python.org/2/library/subprocess.html#module-subprocess)

        Returns
        -------
        Stdout produced after running the program
        """

        try:
            stdout = subprocess.check_output(self.cmd_line, shell=True)
        except subprocess.CalledProcessError as e:
            raise

        return stdout
