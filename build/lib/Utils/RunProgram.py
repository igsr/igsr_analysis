'''
Created on 03 May 2018

@author: ernesto
'''

import os
import subprocess
import pdb
import re

class RunProgram(object):
    '''
    SuperClass used to run an external program within a Python script
    '''


    def __init__(self, program, path=None, args=None, arg_sep=' ', parameters=None,
                 cmd_line=None, downpipe=None, log_name=None, log_file=None):
        '''
        Constructor

        Class variables
        ---------------
        program : str, Required
                  Program to be run
        path : str, Optional
               Folder containing the 'program'
        args : list, Optional
               List of named tuples tuples formed by an argument (or option) and 
               its respective value: i.e. arg=namedtuple('Argument', 'option value')
        parameters : list, Optional
                      List of parameters ['c','d']
        arg_sep : char, Optional
                  char used as a separator between argument and value. 
                  i.e. if '=' then we will get 'a'=1
                  Default is a single whitespace (i.e. ' ')
        cmd_line : str, Optional
                   String with command line to run
        downpipe: list of RunProgram objects, Optional
                  List of RunProgram objects that will be executed in a pipe after
                  self.program has been executed
        log_name: str, Optional
                  Name of the logger
        log_file: str, Optional
                  Path to the file that will be used by the logging library

        '''

        self.path = path
        self.program = program
        self.args = args
        self.arg_sep = arg_sep
        self.parameters = parameters
        self.cmd_line = cmd_line
        self.downpipe = downpipe
        self.log_name = log_name
        self.log_file = log_file

        #create the command line if is None
        if self.cmd_line is None:
            cmd_line=""
            
            if self.path is not None: cmd_line="{0}/".format(self.path)
            
            cmd_line="{0}{1} ".format(cmd_line,self.program)

            # construct the command line
            if self.args is not None:
                for option,parameter in self.args:
                    cmd_line+="{0}{1}{2} ".format(option,self.arg_sep,parameter)

            if self.parameters is not None:
                for param in self.parameters:
                    cmd_line+="{0} ".format(param)

            if downpipe is not None:
                
                for runO in downpipe:
                    cmd_line+="|{0}".format(runO.cmd_line)
        self.cmd_line=cmd_line
        

    def run_popen(self, raise_exc=True):
        '''
        Run self.program using subprocess Popen method
        (see https://docs.python.org/2/library/subprocess.html#module-subprocess)

        Arguments
        ---------
        raise_exc: Bool, optional
                   If true, then raise an Exception when error is found. Default= True

        Returns
        -------
        A tuple containing the STDOUT and STDERR from this program
        '''

        log_f=None

        if self.log_file is not None: log_f=open(self.log_file,'w')

        #execute cmd_line
        p = subprocess.Popen(self.cmd_line, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, bufsize=256, universal_newlines=True)

        #stderr
        patt = re.compile('#* ERROR|Error')
        is_exception=False
        stderr=""
        for line in p.stderr:
            line=str(line.rstrip())
            stderr+=line+"\n"
            if log_f is not None: log_f.write(line+"\n")
            m = patt.match(line)
            if m:
                is_exception=True

        #stdout
        stdout=""
        for line in p.stdout:
            line=str(line.rstrip())
            stdout+=line+"\n"
       
        if is_exception is True and raise_exc is True:
            raise Exception(stderr)


        if log_f is not None:log_f.close()

        return (stdout,stderr,is_exception)

    def run_checkoutput(self):
        '''
        Run self.program using subprocess check_output method
        (see https://docs.python.org/2/library/subprocess.html#module-subprocess)

        Returns
        -------
        Stdout produced after running the program
        '''
        
        try:
            stdout = subprocess.check_output(self.cmd_line, shell=True)
        except subprocess.CalledProcessError as e:
            raise
        
        return stdout

