

import unittest
from unittest import mock

from PyHive.TOPMed.Star import Star
from PyHive.TOPMed.Rsem import Rsem
from PyHive.TOPMed.RnaseqcCounts import RnaseqcCounts
from PyHive.TOPMed.MarkDuplicates import MarkDuplicates
from PyHive.TOPMed.IndexBam import IndexBam
from PyHive.TOPMed.SamtoolsQuickcheck import SamtoolsQuickcheck


class AbstractSingularity():
    def set_singularity_type(self):
        raise NotImplementedError()

    def add_properties(self):
        pass

    def setUp(self):
        self.singularity_type = self.set_singularity_type()
        self.params = {
            'root_output_dir': '/a/working/directory',
            'basename': 'HG001.AAA',
        }
        self.properties = {
        }
        self.add_properties()
        self.params.update({arg: arg for arg in self.singularity_type.CMD_ARGS})
        patches = [
            mock.patch.object(self.singularity_type, "__init__", lambda _self, a, b: None),
            mock.patch.object(self.singularity_type, 'param_is_defined', lambda _self, v: v in self.params),
            mock.patch.object(self.singularity_type, 'param_required', lambda _self, v: self.params[v]),
            mock.patch.object(self.singularity_type, 'param', lambda _self, v: self.params[v]),
        ]
        patches.extend([mock.patch.object(self.singularity_type, k, v) for k, v in self.properties.items()])
        for p in patches:
            self.addCleanup(p.stop)
            p.start()

    def test_commands(self):
        s = self.singularity_type(None, None)
        # PIPELINE
        self.assertIsNotNone(s.PIPELINE)
        # CMD_ARGS
        self.assertIsNotNone(s.get_cmd_args())
        # CMD_KWARGS
        if s.CMD_KWARGS:
            self.assertIn('{ARGS}', s.CMD)
        # CMD
        self.assertIsNotNone(s.make_task_command())
        # FILES
        result = s.get_output_file_list()
        self.assertEqual(len(result), len(s.FILES))


class TestStar(AbstractSingularity, unittest.TestCase):
    def set_singularity_type(self):
        return Star

    def add_properties(self):
        self.params.update({'fastq': ['fastq1', 'fastq2']})


class TestRsem(AbstractSingularity, unittest.TestCase):
    def set_singularity_type(self):

        return Rsem


class TestRnaseqcCounts(AbstractSingularity, unittest.TestCase):
    def set_singularity_type(self):
        return RnaseqcCounts


class TestMarkDuplicates(AbstractSingularity, unittest.TestCase):
    def set_singularity_type(self):
        return MarkDuplicates


class TestIndexBam(AbstractSingularity, unittest.TestCase):
    def set_singularity_type(self):
        return IndexBam


class TestQuickCheck(AbstractSingularity, unittest.TestCase):
    def set_singularity_type(self):
        return SamtoolsQuickcheck

