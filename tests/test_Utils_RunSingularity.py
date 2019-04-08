import unittest
from unittest.mock import patch
import string
import os
import tempfile
import shutil
from Utils.RunSingularity import Singularity, random_suffix

"""
https://medium.com/@george.shuklin/mocking-complicated-init-in-python-6ef9850dd202
"""

class TestSingularity_Setup(unittest.TestCase):
    pass

class TestSingularity_GetOutputDirectory(unittest.TestCase):
    pass


class TestSingularity_GetPrefix(unittest.TestCase):
    def test_prefix_uses_basename(self):
        with patch.object(Singularity, "__init__", lambda x, y, z: None), \
                patch.object(Singularity, 'param_required', lambda _self, v: v):
            s = Singularity(None, None)
            s.PIPELINE = 'pipeline'
            result = s.get_prefix()
        self.assertEqual(result, 'basename.pipeline')

    def test_prefix_does_not_include_colon(self):
        with patch.object(Singularity, "__init__", lambda x, y, z: None), \
                patch.object(Singularity, 'param_required', lambda _self, v: 'root:'+v):
            s = Singularity(None, None)
            s.PIPELINE = 'pipeline'
            result = s.get_prefix()
        self.assertEqual(result, 'root.basename.pipeline')


class TestSingularity_GetOutputFileList(unittest.TestCase):
    pass


class TestSingularity_MakeSingularityCommand(unittest.TestCase):
    def test_error_if_working_dir_not_defined(self):
        with patch.object(Singularity, "__init__", lambda x, y, z: None), \
                patch.object(Singularity, 'param_required', lambda _self, v: v):
            s = Singularity(None, None)
            s.working_dir = None
            with self.assertRaisesRegex(AssertionError, "Working directory has not been defined"):
                s.make_singularity_command()

    def test_singularity_executable_used(self):
        with patch.object(Singularity, "__init__", lambda x, y, z: None), \
                patch.object(Singularity, 'param_required', lambda _self, v: v):
            s = Singularity(None, None)
            s.working_dir = '/a/path/to/working/dir'
            cmd = s.make_singularity_command()
        self.assertTrue(cmd.startswith('singularity_executable'))

    def test_execute_command_in_second_position(self):
        with patch.object(Singularity, "__init__", lambda x, y, z: None), \
                patch.object(Singularity, 'param_required', lambda _self, v: v):
            s = Singularity(None, None)
            s.working_dir = '/a/path/to/working/dir'
            cmd = s.make_singularity_command()
        cmd = cmd.split(' ')
        self.assertEqual(cmd[1], 'exec')

    def test_pwd_in_third_location(self):
        path = '/a/path/to/working/dir'
        with patch.object(Singularity, "__init__", lambda x, y, z: None), \
                patch.object(Singularity, 'param_required', lambda _self, v: v):
            s = Singularity(None, None)
            s.working_dir = path
            cmd = s.make_singularity_command()
        cmd = cmd.split(' ')
        self.assertEqual(cmd[2], '--pwd')
        self.assertEqual(cmd[3], path)

    def test_image_path_in_third_position(self):
        path = '/a/path/to/working/dir'
        with patch.object(Singularity, "__init__", lambda x, y, z: None), \
                patch.object(Singularity, 'param_required', lambda _self, v: v):
            s = Singularity(None, None)
            s.working_dir = path
            cmd = s.make_singularity_command()
        cmd = cmd.split(' ')
        self.assertEqual(cmd[-1], 'singularity_cache/singularity_image')


class TestSingularity_MakeTaskCommand(unittest.TestCase):
    pass

class TestSingularity_GetCmdArgs(unittest.TestCase):
    pass

class TestSingularity_UnpackCmdKwargs(unittest.TestCase):
    pass

class TestSingularity_ExecuteProgram(unittest.TestCase):
    pass

class TestSingularity_OpenWorkingDir(unittest.TestCase):
    pass

class TestSingularity_CloseWorkingDir(unittest.TestCase):
    pass

class TestSingularity_WriteOutput(unittest.TestCase):
    pass


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

#
#
# class TestSingularityUnpackCmdArguments(unittest.TestCase):
#     def setUp(self):
#         self.kwargs = {
#             's_cache': "/foo/bar", 's_image': 'myimage@sha256:a2b4c6d8',
#             'output_directory': '/a/dir', 'prefix': 'a-task'
#         }
#
#     def test_unpack_with_no_values(self):
#         rsing = Singularity(**self.kwargs)
#         self.assertFalse(rsing.CMD_KWARGS)
#         self.assertFalse(rsing.cmd_arguments)
#
#         result = rsing.unpack_cmd_arguments()
#         self.assertEqual(result, '')
#
#     def test_unpack_with_argument_values(self):
#         rsing = Singularity(**self.kwargs)
#         rsing.cmd_arguments={'egg': 'chicken', 'sausage': 'pork'}
#         result = rsing.unpack_cmd_arguments()
#         self.assertEqual(result, '')
#
#     def test_unpack_with_kw_argument_values(self):
#         rsing = Singularity(**self.kwargs)
#         rsing.cmd_arguments={'egg': 'chicken', 'sausage': 'pork'}
#         rsing.CMD_KWARGS = ['egg', 'sausage']
#
#         result = rsing.unpack_cmd_arguments()
#         self.assertEqual(result, '--egg chicken --sausage pork')
#
#
# class TestSingularityOpenWorkingDir(unittest.TestCase):
#     def setUp(self):
#         self.output_dir_obj = tempfile.TemporaryDirectory()
#         self.output_dir = self.output_dir_obj.name
#         self.prefix = 'a-task'
#         self.kwargs = {
#             's_cache': "/foo/bar",
#             's_image': 'myimage@sha256:a2b4c6d8',
#             'prefix': self.prefix
#         }
#
#     def test_working_directory_created(self):
#         rsing = Singularity(output_directory=self.output_dir, **self.kwargs)
#         self.assertIsNone(rsing.working_dir)
#         rsing.open_working_dir()
#         result = rsing.working_dir
#         self.assertIsNotNone(result)
#         self.assertTrue(os.path.exists(result))
#         self.assertEqual(self.output_dir, os.path.dirname(result))
#         self.assertTrue(os.path.basename(result).startswith(self.prefix))
#
#     def test_output_directory_created(self):
#         output = os.path.join(self.output_dir, 'population', 'sample')
#         rsing = Singularity(output_directory=output, **self.kwargs)
#         self.assertIsNone(rsing.working_dir)
#         rsing.open_working_dir()
#         result = rsing.working_dir
#         self.assertIsNotNone(result)
#         self.assertTrue(os.path.exists(result))
#         self.assertEqual(output, os.path.dirname(result))
#         self.assertTrue(os.path.basename(result).startswith(self.prefix))
#
#
# class TestSingularityCloseWorkingDir(unittest.TestCase):
#     def setUp(self):
#         self.output_dir_obj = tempfile.TemporaryDirectory()
#         self.output_dir = self.output_dir_obj.name
#         self.working_dir = os.path.join(self.output_dir, 'working_dir')
#         os.makedirs(self.working_dir)
#
#         self.kwargs = {
#             'output_directory': self.output_dir,
#             's_cache': "/foo/bar",
#             's_image': 'myimage@sha256:a2b4c6d8',
#             'prefix': 'a-task'
#         }
#
#     def test_file_moved_to_output_directory(self):
#         rsing = Singularity(**self.kwargs)
#
#         file_name = 'a_file.txt'
#         working_file = os.path.join(self.working_dir, file_name)
#         output_file = os.path.join(self.output_dir, file_name)
#
#         with open(working_file, 'w') as fh:
#             fh.write('text')
#
#         self.assertTrue(os.path.exists(working_file))
#         self.assertFalse(os.path.exists(output_file))
#
#         rsing.working_dir = self.working_dir
#         result = rsing.close_working_dir()
#
#         self.assertTrue(os.path.exists(output_file))
#         self.assertFalse(os.path.exists(working_file))
#         self.assertEqual(result, [output_file, ])
#         self.assertFalse(os.path.exists(self.working_dir))
#
#     def test_error_when_file_moved_to_output_directory_exists(self):
#         rsing = Singularity(**self.kwargs)
#
#         file_name = 'a_file.txt'
#         working_file = os.path.join(self.working_dir, file_name)
#         output_file = os.path.join(self.output_dir, file_name)
#
#         with open(working_file, 'w') as fh_w:
#             fh_w.write('text')
#         with open(output_file, 'w') as fh_o:
#             fh_o.write('text')
#
#         self.assertTrue(os.path.exists(working_file))
#         self.assertTrue(os.path.exists(output_file))
#
#         rsing.working_dir = self.working_dir
#
#         with self.assertRaisesRegex(shutil.Error,
#                                     "Cannot close working directory, files already exist"):
#             rsing.close_working_dir()
#
#
