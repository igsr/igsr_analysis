import unittest
import os
import tempfile
import shutil
from Utils.RunSingularity import Singularity


class TestSingularity(unittest.TestCase):
    def test(self):
        pass


class TestSingularityMakeSingularityCmd(unittest.TestCase):
    def setUp(self):
        self.kwargs = {
            'output_directory': '/a/dir', 'prefix': 'a-task'
        }

    def test_default_singularity_executable_used(self):
        rsing = Singularity(s_cache='', s_image='', **self.kwargs)
        cmd = rsing.make_singularity_command()
        self.assertTrue(cmd.startswith('singularity'))

    def test_specified_singularity_executable_used(self):
        exe = "/bin/singularity"
        rsing = Singularity(s_cache='', s_image='', s_executable=exe, **self.kwargs)
        cmd = rsing.make_singularity_command()
        self.assertTrue(cmd.startswith(exe))

    def test_execute_command_in_second_position(self):
        rsing = Singularity(s_cache='', s_image='myimage@sha256:a2b4c6d8', **self.kwargs)
        cmd = rsing.make_singularity_command()
        cmd = cmd.split(' ')
        self.assertEqual(cmd[1], 'exec')

    def test_image_path_in_third_position(self):
        cache = "/foo/bar"
        img = 'myimage@sha256:a2b4c6d8'
        rsing = Singularity(s_cache=cache, s_image=img, **self.kwargs)
        cmd = rsing.make_singularity_command()
        cmd = cmd.split(' ')
        self.assertEqual(cmd[2], cache+'/'+img)

    def test_working_dir_in_command(self):
        cache = "/foo/bar"
        img = 'myimage@sha256:a2b4c6d8'
        rsing = Singularity(s_cache=cache, s_image=img, **self.kwargs)
        rsing.working_dir = os.path.join(self.kwargs['output_directory'], 'working_a1b2c')
        self.assertIn(f"--pwd {rsing.working_dir}", rsing.make_singularity_command())

    def test_error_if_working_dir_not_defined(self):
        cache = "/foo/bar"
        img = 'myimage@sha256:a2b4c6d8'
        rsing = Singularity(s_cache=cache, s_image=img, **self.kwargs)
        with self.assertRaisesRegex(AssertionError, "Working directory has not been defined"):
            rsing.make_singularity_command()



class TestSingularityUnpackCmdArguments(unittest.TestCase):
    def setUp(self):
        self.kwargs = {
            's_cache': "/foo/bar", 's_image': 'myimage@sha256:a2b4c6d8',
            'output_directory': '/a/dir', 'prefix': 'a-task'
        }

    def test_unpack_with_no_values(self):
        rsing = Singularity(**self.kwargs)
        self.assertFalse(rsing.CMD_KWARGS)
        self.assertFalse(rsing.cmd_arguments)

        result = rsing.unpack_cmd_arguments()
        self.assertEqual(result, '')

    def test_unpack_with_argument_values(self):
        rsing = Singularity(**self.kwargs)
        rsing.cmd_arguments={'egg': 'chicken', 'sausage': 'pork'}
        result = rsing.unpack_cmd_arguments()
        self.assertEqual(result, '')

    def test_unpack_with_kw_argument_values(self):
        rsing = Singularity(**self.kwargs)
        rsing.cmd_arguments={'egg': 'chicken', 'sausage': 'pork'}
        rsing.CMD_KWARGS = ['egg', 'sausage']

        result = rsing.unpack_cmd_arguments()
        self.assertEqual(result, '--egg chicken --sausage pork')


class TestSingularityOpenWorkingDir(unittest.TestCase):
    def setUp(self):
        self.output_dir_obj = tempfile.TemporaryDirectory()
        self.output_dir = self.output_dir_obj.name
        self.prefix = 'a-task'
        self.kwargs = {
            's_cache': "/foo/bar",
            's_image': 'myimage@sha256:a2b4c6d8',
            'prefix': self.prefix
        }

    def test_working_directory_created(self):
        rsing = Singularity(output_directory=self.output_dir, **self.kwargs)
        self.assertIsNone(rsing.working_dir)
        rsing.open_working_dir()
        result = rsing.working_dir
        self.assertIsNotNone(result)
        self.assertTrue(os.path.exists(result))
        self.assertEqual(self.output_dir, os.path.dirname(result))
        self.assertTrue(os.path.basename(result).startswith(self.prefix))

    def test_output_directory_created(self):
        output = os.path.join(self.output_dir, 'population', 'sample')
        rsing = Singularity(output_directory=output, **self.kwargs)
        self.assertIsNone(rsing.working_dir)
        rsing.open_working_dir()
        result = rsing.working_dir
        self.assertIsNotNone(result)
        self.assertTrue(os.path.exists(result))
        self.assertEqual(output, os.path.dirname(result))
        self.assertTrue(os.path.basename(result).startswith(self.prefix))


class TestSingularityCloseWorkingDir(unittest.TestCase):
    def setUp(self):
        self.output_dir_obj = tempfile.TemporaryDirectory()
        self.output_dir = self.output_dir_obj.name
        self.working_dir = os.path.join(self.output_dir, 'working_dir')
        os.makedirs(self.working_dir)

        self.kwargs = {
            'output_directory': self.output_dir,
            's_cache': "/foo/bar",
            's_image': 'myimage@sha256:a2b4c6d8',
            'prefix': 'a-task'
        }

    def test_file_moved_to_output_directory(self):
        rsing = Singularity(**self.kwargs)

        file_name = 'a_file.txt'
        working_file = os.path.join(self.working_dir, file_name)
        output_file = os.path.join(self.output_dir, file_name)

        with open(working_file, 'w') as fh:
            fh.write('text')

        self.assertTrue(os.path.exists(working_file))
        self.assertFalse(os.path.exists(output_file))

        rsing.working_dir = self.working_dir
        result = rsing.close_working_dir()

        self.assertTrue(os.path.exists(output_file))
        self.assertFalse(os.path.exists(working_file))
        self.assertEqual(result, [output_file, ])
        self.assertFalse(os.path.exists(self.working_dir))

    def test_error_when_file_moved_to_output_directory_exists(self):
        rsing = Singularity(**self.kwargs)

        file_name = 'a_file.txt'
        working_file = os.path.join(self.working_dir, file_name)
        output_file = os.path.join(self.output_dir, file_name)

        with open(working_file, 'w') as fh_w:
            fh_w.write('text')
        with open(output_file, 'w') as fh_o:
            fh_o.write('text')

        self.assertTrue(os.path.exists(working_file))
        self.assertTrue(os.path.exists(output_file))

        rsing.working_dir = self.working_dir

        with self.assertRaisesRegex(shutil.Error,
                                    "Cannot close working directory, files already exist"):
            rsing.close_working_dir()


