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


    def __init__(self, program, path=None, args=None, parameters=None,
                 cmd_line=None, downpipe=None):
        '''
        Constructor

        Class variables
        ---------------
        program : str, Required
                  Program to be run
        path : str, Optional
               Folder containing the 'program'
        args : list, Optional
               List of tuples formed by an argument (or option) and 
               its respective value: i.e. [(-a, 20),(-b, 'aname')]
        parameters : list, Optional
                      List of parameters ['c','d']
        cmd_line : str, Optional
                   String with command line to run
        downpipe: list of RunProgram objects, Optional
                  List of RunProgram objects that will be executed in a pipe after
                  self.program has been executed

        '''

        self.path = path
        self.program = program
        self.args = args
        self.parameters = parameters
        self.cmd_line = cmd_line
        self.downpipe = downpipe

        #create the command line if is None
        if self.cmd_line is None:
            cmd_line=""
            
            if self.path is not None: cmd_line="{0}/".format(self.path)
            
            cmd_line="{0}{1}".format(cmd_line,self.program)

            # construct the command line
            if self.args is not None:
                for option,parameter in self.args:
                    cmd_line="{0} {1} {2} ".format(cmd_line,option,parameter)

            if self.parameters is not None:
                for param in self.parameters:
                    cmd_line="{0} {1} ".format(cmd_line,param)
                    
            if downpipe is not None:
                for runO in downpipe:
                    cmd_line="{0} | {1}".format(cmd_line,runO.cmd_line)
        self.cmd_line=cmd_line
        

    def run(self):
        '''
        Run self.program

        Returns
        -------
        A tuple containing the STDOUT and STDERR from this program
        '''

        #execute cmd_line
        try:
            p = subprocess.Popen(self.cmd_line, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,universal_newlines=True)
            stdout, stderr = p.communicate()
        except subprocess.CalledProcessError as exp:
            print("Something went wrong while {0}".format(self.program))
            raise Exception(exp.output)

        return (stdout,stderr)

            
