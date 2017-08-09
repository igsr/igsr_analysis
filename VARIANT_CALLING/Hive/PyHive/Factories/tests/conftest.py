import pytest

def pytest_addoption(parser):
    parser.addoption('--faix', default='data/canonical_chros.fa.fai' ,action='store_true', help='Path to faix file')
    parser.addoption('--hive_lib', default='~/lib/ensembl-hive_2.4/' ,action='store_true', help='Path folder containing eHive scripts')

