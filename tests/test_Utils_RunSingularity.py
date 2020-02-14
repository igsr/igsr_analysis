import unittest
from unittest import mock
import string
import os
import re
import tempfile
import subprocess
import shutil
from Utils.RunSingularity import Singularity, random_suffix

"""
https://medium.com/@george.shuklin/mocking-complicated-init-in-python-6ef9850dd202
"""

class TestSingularity_Setup(unittest.TestCase):
    pass

class TestSingularity_SessionSuffix(unittest.TestCase):
    def setUp(self):
        patches = [
            mock.patch.object(Singularity, "__init__", lambda _self, a, b: None)
        ]
        for p in patches:
            self.addCleanup(p.stop)
            p.start()

    def test_session_suffix_is_randomised(self):
        results = []
        for __ in range(20):
            s = Singularity(None, None)
            results.append(s.session_suffix)
        self.assertEqual(len(results), len(set(results)))

    def test_session_suffix_can_be_recalled(self):
        s = Singularity(None, None)
        result = {s.session_suffix for __ in range(20)}
        self.assertEqual(len(result), 1)


class TestSingularity_GeneralProperties(unittest.TestCase):
    def setUp(self):
        self.params = {
            'basename': 'base:name',
        }
        patches = [
            mock.patch.object(Singularity, "__init__", lambda _self, a, b: None),
            mock.patch.object(Singularity, 'param_exists', lambda _self, v: v in self.params),
            mock.patch.object(Singularity, 'param_required', lambda _self, v: self.params[v]),
            mock.patch.object(Singularity, 'param', lambda _self, v: self.params[v]),
        ]
        self.properties = {
            'PIPELINE': 'mock_pipeline',
            'session_suffix': 'abcd123',
        }
        patches.extend([mock.patch.object(Singularity, k, v) for k, v in self.properties.items()])
        for p in patches:
            self.addCleanup(p.stop)
            p.start()

    # ==================================================================== #
    # Base Prefix
    def test_base_prefix_uses_basename(self):
        self.params['basename'] = 'a_different_basename'
        s = Singularity(None, None)
        result = s.base_prefix
        self.assertEqual(result, self.params['basename'])

    def test_base_prefix_replaces_colon(self):
        expected = self.params['basename'] .replace(':', '.')
        self.assertNotEqual(expected, self.params['basename'])
        s = Singularity(None, None)
        result = s.base_prefix
        self.assertEqual(result, expected)

    # ==================================================================== #
    # Prefix
    def test_prefix_uses_base_prefix(self):
        s = Singularity(None, None)
        result = s.prefix
        self.assertTrue(result.startswith(s.base_prefix))

    def test_prefix_uses_pipeline_name(self):
        s = Singularity(None, None)
        result = s.prefix
        self.assertTrue(result.endswith(self.properties['PIPELINE']))


class TestSingularity_DirectoryProperties(unittest.TestCase):
    def setUp(self):
        self.params = {
            'root_output_dir': '/root/output/dir',
            'basename': 'base:name',
        }
        patches = [
            mock.patch.object(Singularity, "__init__", lambda _self, a, b: None),
            mock.patch.object(Singularity, 'param_is_defined', lambda _self, v: v in self.params),
            mock.patch.object(Singularity, 'param_required', lambda _self, v: self.params[v]),
            mock.patch.object(Singularity, 'param', lambda _self, v: self.params[v]),
        ]
        self.properties = {
            'PIPELINE': 'mock_pipeline',
            'session_suffix': 'abcd123',
        }
        patches.extend([mock.patch.object(Singularity, k, v) for k, v in self.properties.items()])
        for p in patches:
            self.addCleanup(p.stop)
            p.start()

    # ==================================================================== #
    # Output Directory
    def test_output_is_root_if_no_dir_labels(self):
        s = Singularity(None, None)
        result = s.output_directory
        self.assertEqual(result, self.params['root_output_dir'])

    def test_output_is_nested_if_dir_labels(self):
        self.params.update({
            'dir_label_params': ['a', 'b', 'c'],
            'a': 'a_dir',
            'b': 'b_dir',
            'c': 'c_dir'
        })
        expected = os.path.join(*[self.params[k] for k in
                                  ('root_output_dir', 'a', 'b', 'c')])
        s = Singularity(None, None)
        result = s.output_directory
        self.assertEqual(result, expected)

    # ==================================================================== #
    # Working Directory
    def test_working_directory(self):
        output_dir = '/a/root/outptut/dir/for/sample'
        prefix = 'a.prefix'
        expected = os.path.join(output_dir, "{}_{}".format(prefix, self.properties['session_suffix']))
        with mock.patch.object(Singularity, 'output_directory', output_dir),\
                mock.patch.object(Singularity, 'prefix', prefix):
            s = Singularity(None, None)
            result = s.working_directory
        self.assertEqual(result, expected)

    # ==================================================================== #
    # Logging Directory
    def test_log_directory(self):
        output_dir = '/a/root/outptut/dir/for/sample'
        expected = os.path.join(output_dir, 'ehive_log')
        with mock.patch.object(Singularity, 'output_directory', output_dir):
            s = Singularity(None, None)
            result = s.log_directory
        self.assertEqual(result, expected)

    # ==================================================================== #
    # Logging Basename
    def test_log_basename_directory(self):
        output_dir = '/a/root/outptut/dir/for/sample'
        expected = os.path.join(output_dir, 'ehive_log', 'base.name.mock_pipeline.abcd123')
        with mock.patch.object(Singularity, 'output_directory', output_dir):
            s = Singularity(None, None)
            result = s.log_basename
        self.assertEqual(result, expected)


class TestSingularity_GetOutputFileList(unittest.TestCase):
    def setUp(self):
        self.properties = {
            'output_directory': '/an/output/directory',
            'prefix': 'HG001.AAA.a_pipeline',
        }
        patches = [mock.patch.object(Singularity, "__init__", lambda _self, a, b: None)]
        patches.extend([mock.patch.object(Singularity, k, v) for k, v in self.properties.items()])
        for p in patches:
            self.addCleanup(p.stop)
            p.start()

    def test_files_without_replacement(self):
        files = {'foo': 'bar'}
        expected = {k: os.path.join(self.properties['output_directory'], v) for k, v in files.items()}
        with mock.patch.object(Singularity, 'FILES', files):
            s = Singularity(None, None)
            result = s.get_output_file_list()
        self.assertEqual(result, expected)

    def test_files_with_prefix_added(self):
        files = {'foo': '{PREFIX}.bar'}
        expected = {k: os.path.join(self.properties['output_directory'], v.format(PREFIX=self.properties['prefix']))
                    for k, v in files.items()}
        with mock.patch.object(Singularity, 'FILES', files):
            s = Singularity(None, None)
            result = s.get_output_file_list()
        self.assertEqual(result, expected)

    def test_error_on_files_with_non_prefix_format_string(self):
        files = {'foo': '{bar}.bar'}
        with mock.patch.object(Singularity, 'FILES', files):
            s = Singularity(None, None)
            with self.assertRaises(KeyError):
                result = s.get_output_file_list()


class TestSingularity_MakeSingularityCommand(unittest.TestCase):
    def setUp(self):
        self.working_dir = '/a/path/to/working/dir'
        patches = [
            mock.patch.object(Singularity, "__init__", lambda _self, a, b: None),
            mock.patch.object(Singularity, 'param_required', lambda _self, v: v),
            mock.patch.object(Singularity, 'working_directory', self.working_dir)
        ]
        for p in patches:
            self.addCleanup(p.stop)
            p.start()

    def test_singularity_executable_used(self):
        s = Singularity(None, None)
        cmd = s.make_singularity_command()
        self.assertTrue(cmd.startswith('singularity_executable'))

    def test_execute_command_in_second_position(self):
        s = Singularity(None, None)
        cmd = s.make_singularity_command()
        cmd = cmd.split(' ')
        self.assertEqual(cmd[1], 'exec')

    def test_pwd_in_third_location(self):
        s = Singularity(None, None)
        cmd = s.make_singularity_command()
        cmd = cmd.split(' ')
        self.assertEqual(cmd[2], '--pwd')
        self.assertEqual(cmd[3], self.working_dir)

    def test_image_path_in_third_position(self):
        s = Singularity(None, None)
        cmd = s.make_singularity_command()
        cmd = cmd.split(' ')
        self.assertEqual(cmd[-1], 'singularity_cache/singularity_image')


class TestSingularity_MakeTaskCommand(unittest.TestCase):
    def setUp(self):
        self.params = {}
        patches = [
            mock.patch.object(Singularity, "__init__", lambda _self, a, b: None),
            mock.patch.object(Singularity, 'param_is_defined', lambda _self, v: v in self.params),
            mock.patch.object(Singularity, 'param_required', lambda _self, v: self.params[v]),
            mock.patch.object(Singularity, 'param', lambda _self, v: self.params[v]),
        ]
        self.properties = {
            'working_directory': '/a/working/directory',
            'prefix': 'HG001.AAA.a_pipeline',
        }
        patches.extend([mock.patch.object(Singularity, k, v) for k, v in self.properties.items()])
        for p in patches:
            self.addCleanup(p.stop)
            p.start()

    def test_cmd_with_no_replacements(self):
        cmd = "echo $PATH"
        with mock.patch.object(Singularity, 'CMD', cmd):
            s = Singularity(None, None)
            result = s.make_task_command()
        self.assertEqual(result, cmd)

    def test_cmd_with_working_directory_replacement(self):
        cmd = "ls {WORKING_DIR}"
        expected = cmd.format(WORKING_DIR=self.properties['working_directory'])
        with mock.patch.object(Singularity, 'CMD', cmd):
            s = Singularity(None, None)
            result = s.make_task_command()
        self.assertEqual(result, expected)

    def test_cmd_with_unfulfilled_args_option(self):
        cmd_kwargs = ["optiona", 'optionb']
        cmd = "command {ARGS}"
        expected = cmd.format(ARGS='')
        with mock.patch.object(Singularity, 'CMD', cmd), \
                mock.patch.object(Singularity, 'CMD_KWARGS', cmd_kwargs):
            s = Singularity(None, None)
            result = s.make_task_command()
        self.assertEqual(result, expected)

    def test_cmd_with_fulfilled_args_option(self):
        cmd_kwargs = ["optiona", 'optionb']
        self.params['optionb'] = 'valueb'
        cmd = "command {ARGS}"
        expected = cmd.format(ARGS='--optionb valueb')
        with mock.patch.object(Singularity, 'CMD', cmd), \
                mock.patch.object(Singularity, 'CMD_KWARGS', cmd_kwargs):
            s = Singularity(None, None)
            result = s.make_task_command()
        self.assertEqual(result, expected)

    def test_cmd_with_multiple_fulfilled_args_option(self):
        cmd_kwargs = ["optiona", 'optionb']
        self.params['optiona'] = 'valuea'
        self.params['optionb'] = 'valueb'
        cmd = "command {ARGS}"
        expected = cmd.format(ARGS='--optiona valuea --optionb valueb')
        with mock.patch.object(Singularity, 'CMD', cmd), \
                mock.patch.object(Singularity, 'CMD_KWARGS', cmd_kwargs):
            s = Singularity(None, None)
            result = s.make_task_command()
        self.assertEqual(result, expected)


class TestSingularity_GetCmdArgs(unittest.TestCase):
    def setUp(self):
        self.params = {}
        patches = [
            mock.patch.object(Singularity, "__init__", lambda _self, a, b: None),
            mock.patch.object(Singularity, 'param_is_defined', lambda _self, v: v in self.params),
            mock.patch.object(Singularity, 'param_required', lambda _self, v: self.params[v]),
            mock.patch.object(Singularity, 'param', lambda _self, v: self.params[v]),
        ]
        self.properties = {
            'working_directory': '/a/working/directory',
            'prefix': 'HG001.AAA.a_pipeline',
        }
        patches.extend([mock.patch.object(Singularity, k, v) for k, v in self.properties.items()])
        for p in patches:
            self.addCleanup(p.stop)
            p.start()


class TestSingularity_UnpackCmdKwargs(unittest.TestCase):
    def setUp(self):
        self.params = {}
        patches = [
            mock.patch.object(Singularity, "__init__", lambda _self, a, b: None),
            mock.patch.object(Singularity, 'param_is_defined', lambda _self, v: v in self.params),
            mock.patch.object(Singularity, 'param_required', lambda _self, v: self.params[v]),
            mock.patch.object(Singularity, 'param', lambda _self, v: self.params[v]),
        ]
        for p in patches:
            self.addCleanup(p.stop)
            p.start()

    def test_unpack_with_no_kwargs(self):
        cmd_kwargs = []
        self.params['egg'] = 'scrambled'
        expected = ''
        with mock.patch.object(Singularity, 'CMD_KWARGS', cmd_kwargs):
            s = Singularity(None, None)
            result = s.unpack_cmd_kwargs()
        self.assertEqual(result, expected)

    def test_unpack_with_no_matches(self):
        cmd_kwargs = ["egg", 'sausage']
        expected = ''
        with mock.patch.object(Singularity, 'CMD_KWARGS', cmd_kwargs):
            s = Singularity(None, None)
            result = s.unpack_cmd_kwargs()
        self.assertEqual(result, expected)

    def test_cmd_with_multiple_fulfilled_args_option(self):
        cmd_kwargs = ["egg", 'sausage']
        self.params['egg'] = 'poached'
        self.params['sausage'] = 'fried'
        expected = '--egg poached --sausage fried'
        with mock.patch.object(Singularity, 'CMD_KWARGS', cmd_kwargs):
            s = Singularity(None, None)
            result = s.unpack_cmd_kwargs()
        self.assertEqual(result, expected)


class TestSingularity_ExecuteProgram(unittest.TestCase):
    def setUp(self):
        patches = [
            mock.patch.object(Singularity, "__init__", lambda _self, a, b: None),
        ]
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)

        self.log_base = os.path.join(self.temp_dir.name, 'LOGGING')
        self.properties = {'log_basename': self.log_base}
        patches.extend([mock.patch.object(Singularity, k, v) for k, v in self.properties.items()])
        for p in patches:
            self.addCleanup(p.stop)
            p.start()

    def test_log_stored_for_executed_program(self):
        cmd = 'expr 100 + 10'
        s = Singularity(None, None)
        s.execute_program(cmd)

        self.assertTrue(os.path.exists(f'{self.log_base}.stderr.txt'))
        with open(f'{self.log_base}.stderr.txt', 'r') as fh:
            self.assertFalse(fh.read())

        self.assertTrue(os.path.exists(f'{self.log_base}.stdout.txt'))
        with open(f'{self.log_base}.stdout.txt', 'r') as fh:
            self.assertEqual(fh.read(), '110\n')

    def test_error_log_stored_for_executed_program(self):
        cmd = 'expr 100 + 10; foobar'
        s = Singularity(None, None)

        err_msg = f"Command '{cmd}' returned non-zero exit status 127."

        with self.assertRaisesRegex(subprocess.CalledProcessError,
                                    re.escape(err_msg)):
            s.execute_program(cmd)

        self.assertTrue(os.path.exists(f'{self.log_base}.stderr.txt'))
        with open(f'{self.log_base}.stderr.txt', 'r') as fh:
            self.assertEqual(fh.read(), f"/bin/sh: foobar: command not found\n\nERROR: {err_msg}\n")

        self.assertTrue(os.path.exists(f'{self.log_base}.stdout.txt'))
        with open(f'{self.log_base}.stdout.txt', 'r') as fh:
            self.assertEqual(fh.read(), f"110\n\nERROR: {err_msg}\n")


class TestSingularity_CloseWorkingDir(unittest.TestCase):
    def setUp(self):
        patches = [
            mock.patch.object(Singularity, "__init__", lambda _self, a, b: None),
        ]

        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)

        self.properties = {'output_directory': os.path.join(self.temp_dir.name, 'output_dir')}
        self.properties['working_directory'] = os.path.join(self.properties['output_directory'], 'temporary_dir')
        for v in self.properties.values():
            os.makedirs(v, exist_ok=True)
        patches.extend([mock.patch.object(Singularity, k, v) for k, v in self.properties.items()])
        for p in patches:
            self.addCleanup(p.stop)
            p.start()

    def test_working_dir_with_no_content(self):
        self.assertTrue(os.path.exists(self.properties['output_directory']))
        self.assertTrue(os.path.exists(self.properties['working_directory']))

        s = Singularity(None, None)
        s.close_working_dir()

        self.assertTrue(os.path.exists(self.properties['output_directory']))
        self.assertFalse(os.path.exists(self.properties['working_directory']))

    def test_working_dir_with_content(self):
        test_file = 'test.txt'
        self.assertTrue(os.path.exists(self.properties['output_directory']))
        self.assertTrue(os.path.exists(self.properties['working_directory']))

        with open(os.path.join(self.properties['working_directory'], test_file), 'w') as fh:
            fh.write(test_file)

        s = Singularity(None, None)
        s.close_working_dir()

        self.assertTrue(os.path.exists(self.properties['output_directory']))
        self.assertFalse(os.path.exists(self.properties['working_directory']))
        self.assertTrue(os.path.exists(os.path.join(self.properties['output_directory'], test_file)))

    def test_error_if_file_already_exists(self):
        test_file = 'test.txt'
        self.assertTrue(os.path.exists(self.properties['output_directory']))
        self.assertTrue(os.path.exists(self.properties['working_directory']))

        with open(os.path.join(self.properties['working_directory'], test_file), 'w') as fh:
            fh.write(test_file)
        self.assertTrue(os.path.exists(os.path.join(self.properties['working_directory'], test_file)))

        with open(os.path.join(self.properties['output_directory'], test_file), 'w') as fh:
            fh.write(test_file)
        self.assertTrue(os.path.exists(os.path.join(self.properties['output_directory'], test_file)))

        s = Singularity(None, None)
        with self.assertRaises(shutil.Error):
            s.close_working_dir()


class TestSingularity_WriteOutput(unittest.TestCase):
    def setUp(self):
        patches = [
            mock.patch.object(Singularity, "__init__", lambda _self, a, b: None),
        ]

        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)

        self.properties = {
            'output_directory': os.path.join(self.temp_dir.name, 'output_dir'),
            'prefix': 'HG0001.ABC'
        }
        for v in self.properties.values():
            os.makedirs(v, exist_ok=True)
        patches.extend([mock.patch.object(Singularity, k, v) for k, v in self.properties.items()])
        for p in patches:
            self.addCleanup(p.stop)
            p.start()

    def test_write_output_with_no_files(self):
        with mock.patch.object(Singularity, "dataflow") as patch_dataflow, \
                mock.patch.object(Singularity, "warning") as patch_warning:
            s = Singularity(None, None)
            s.write_output()
        patch_dataflow.assert_called_once_with(dict(), 1)

    def test_write_output_with_content(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            files = {'file_a': os.path.join(temp_dir, 'file_a.txt'),
                     'file_b': os.path.join(temp_dir, 'file_b.txt')}
            for name, path in files.items():
                with open(path, 'w') as fh:
                    fh.write(name)

            with mock.patch.object(Singularity, "dataflow") as patch_dataflow, \
                    mock.patch.object(Singularity, "warning") as patch_warning, \
                    mock.patch.object(Singularity, "get_output_file_list", lambda *args, **kwargs: files) as patch_output_file_list:
                s = Singularity(None, None)
                s.write_output()
            patch_dataflow.assert_called_once_with(files, 1)

    def test_error_with_missing_file(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            files = {'file_a': os.path.join(temp_dir, 'file_a.txt'),
                     'file_b': os.path.join(temp_dir, 'file_b.txt')}

            with open(files['file_a'], 'w') as fh:
                fh.write('file_a')
            missing_file = files['file_b']
            self.assertFalse(os.path.exists(missing_file))

            with mock.patch.object(Singularity, "dataflow") as patch_dataflow, \
                    mock.patch.object(Singularity, "warning") as patch_warning, \
                    mock.patch.object(Singularity, "get_output_file_list", lambda *args, **kwargs: files) as patch_output_file_list:
                s = Singularity(None, None)
                with self.assertRaisesRegex(AssertionError, re.escape(f"Missing file: {missing_file}")):
                    s.write_output()


class TestRandomSuffix(unittest.TestCase):
    def test_suffix_length(self):
        res = random_suffix()
        self.assertEqual(len(res), 5)

        l = 3
        res = random_suffix(l)
        self.assertEqual(len(res), l)

        l = 8
        res = random_suffix(l)
        self.assertEqual(len(res), l)

    def test_suffix_lowercase_and_numerical(self):
        valid = set(string.ascii_lowercase + string.digits)
        for i in range(50):
            suffix = random_suffix(20)
            self.assertTrue(all(s in valid for s in suffix))

