import os
import pytest
import pdb

from Utils.RunProgram import RunProgram 

# test_Utils_RunProgram.py


def test_run_program_withoneparam():
    '''
    Test function for the 'run' method

    It will run the shell 'echo' command with one parameter
    '''
    
    runner=RunProgram(program='echo',parameters=['hello'])

    stdout,stderr=runner.run()

    assert stdout=="hello\n"

def test_run_program_withstderr():
    '''
    Test function for the 'run' method

    It will try to run the shell 'echo' incorrectly and will dump some content on STDERR 
    '''

    runner=RunProgram(program='cho', parameters=['-n','hello'])

    stdout,stderr=runner.run()

    assert stderr is not ""

def test_run_inapipe():
    '''
    Test function for the 'run' method with some additional commands executed in a pipe
    '''

    down_cmd=RunProgram(program='wc', parameters=['-l'])
    
    runner=RunProgram(program='cat', parameters=['data/newheader.txt'],downpipe=[down_cmd])

    stdout,stderr=runner.run()

    assert stdout=='67\n'
