import unittest
import pytest
import pdb

from Utils import RunProgram

# test_Utils_RunProgram.py

def test_run_program_withoneparam():
    """
    Test function for the 'run' method

    It will run the shell 'echo' command with one parameter
    """

    runner = RunProgram(program='echo', parameters=['hello'])

    stdout = runner.run_checkoutput()

    assert stdout.decode("utf-8") == "hello\n"

def test_run_program_withstderr():
    """
    Test function for the 'run' method

    It will try to run the shell 'echo' incorrectly and will throw an Exception
    """

    runner = RunProgram(program='cho', parameters=['-n', 'hello'])

    with pytest.raises(Exception):
        runner.run_checkoutput()

def test_run_inapipe():
    """
    Test function for the 'run' method with some additional commands executed in a pipe
    """

    down_cmd = RunProgram(program='wc', parameters=['-l'])

    runner = RunProgram(program='cat', parameters=['data/newheader.txt'], downpipe=[down_cmd])

    stdout = runner.run_checkoutput()

    assert stdout.decode("utf-8").replace(' ', '') == '67\n'

def test_raises_Excp1():
    """
    This test returns exception because of conflicting instructions on how to use
    the dependency. self.path and self.use_docker=True are mutually conflict
    """
    with pytest.raises(Exception) as e_info:
        runner = RunProgram(program='cat', settingf='./data/settings.ini',
                            use_docker=True, path='~/bin')

def test_run_program_w_docker():
    """
    Test to run a containerized application (BCFTools in this case)
    """
    runner = RunProgram(program='bcftools',
                        parameters=['--version-only'],
                        settingf='./data/settings.ini',
                        use_docker=True)

    stdout = runner.run_checkoutput()
    stdout_str = stdout.decode("utf-8").rstrip("\n")
    # check output of running bcftools --version-only
    assert '1.9+htslib-1.9' == '1.9+htslib-1.9'

def test_run_program_w_docker1():
    """
    Test to run a containerized application (BCFTools in this case)
    with an argument
    """
    runner = RunProgram(program='bcftools',
                        parameters=['view', '-h', 'data/GLs.HG00136.correct.chr22:20000085-20010512.beagle.vcf.gz'],
                        settingf='./data/settings.ini',
                        use_docker=True)

    stdout = runner.run_checkoutput()
    stdout_str = stdout.decode("utf-8").split("\n")[0]
    assert '##fileformat=VCFv4.2' == stdout_str

def test_run_program_w_docker2():
    """
    Test to run 2 containerized applications (BCFTools in this case)
    in a pipe
    """
    runner2 = RunProgram(program='bcftools',
                         parameters=['view', '-H'],
                         settingf='./data/settings.ini')

    pipelist = [runner2]

    runner1 = RunProgram(program='bcftools',
                         parameters=['view', 'data/GLs.HG00136.correct.chr22:20000085-20010512.beagle.vcf.gz'],
                         settingf='./data/settings.ini',
                         downpipe=pipelist,
                         use_docker=True)

    stdout = runner1.run_checkoutput()
    stdout_str = stdout.decode("utf-8").split("\n")[0]
    assert '22\t20000086\trs138720731\tT\tC\t.\tPASS\tAR2=0;DR2=0;AF=0\tGT:DS:GP\t0/0:0:1,0,0' == stdout_str

class TestRunProgramInit(unittest.TestCase):
    def test_error_on_no_cmdline_or_program(self):
        with self.assertRaisesRegex(ValueError, "Parameter 'cmd_line'"
                                                " or 'program' must be provided."):
            RunProgram(program=None, cmd_line=None)

class TestRunProgramCreateCommandLine(unittest.TestCase):
    def test_command_line_formation_w_program(self):
        kwargs = {'program': 'program.exe'}
        rp = RunProgram(**kwargs)
        cmd_result = rp.create_command_line()
        self.assertEqual(cmd_result, "{program}".format(**kwargs))

    def test_command_line_formation_w_path_program(self):
        kwargs = {'path': '/a/path', 'program': 'program.exe'}
        rp = RunProgram(**kwargs)
        cmd_result = rp.create_command_line()
        self.assertEqual(cmd_result, "{path}/{program}".format(**kwargs))

    def test_command_line_formation_w_path_program_trailing_slash(self):
        kwargs = {'path': '/a/path/', 'program': 'program.exe'}
        rp = RunProgram(**kwargs)
        cmd_result = rp.create_command_line()
        self.assertEqual(cmd_result, "{path}{program}".format(**kwargs))

    def test_command_line_formation_w_args(self):
        ak, av, bk, bv = 'chicken', 'egg', 'pork', 'sausage'
        kwargs = {'program': 'program.exe', 'args': [(ak, av), (bk, bv)]}
        rp = RunProgram(**kwargs)
        cmd_result = rp.create_command_line()
        self.assertEqual(cmd_result, f"{kwargs['program']} {ak} {av} {bk} {bv}")

    def test_command_line_formation_w_args_sep(self):
        ak, av, bk, bv = 'chicken', 'egg', 'pork', 'sausage'
        kwargs = {
            'program': 'program.exe',
            'args': [(ak, av), (bk, bv)],
            'arg_sep': '='
        }
        rp = RunProgram(**kwargs)
        cmd_result = rp.create_command_line()
        self.assertEqual(cmd_result, f"{kwargs['program']} {ak}={av} {bk}={bv}")

    def test_command_line_formation_w_parameters(self):
        av, bv = 'egg', 'sausage'
        kwargs = {'program': 'program.exe', 'parameters': [av, bv]}
        rp = RunProgram(**kwargs)
        cmd_result = rp.create_command_line()
        self.assertEqual(cmd_result, f"{kwargs['program']} {av} {bv}")

    def test_command_line_formation_w_downpipe(self):
        a, b = 'egg', 'sausage'
        a_rp, b_rp = RunProgram(cmd_line=a), RunProgram(cmd_line=b)
        kwargs = {'program': 'program.exe', 'downpipe': [a_rp, b_rp]}
        rp = RunProgram(**kwargs)
        cmd_result = rp.create_command_line()
        self.assertEqual(cmd_result, f"{kwargs['program']} | {a} | {b}")