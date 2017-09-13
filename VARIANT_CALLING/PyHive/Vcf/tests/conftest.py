import pytest

def pytest_addoption(parser):
    parser.addoption('--vcf', default='data/test.vcf.gz' ,action='store_true', help='Path to vcf file')
    parser.addoption('--bcftools_folder', default='~/bin/bcftools/' ,action='store_true', help='Folder containing bcftools binaries')
    parser.addoption('--hive_lib', default='~/lib/ensembl-hive_2.4/' ,action='store_true', help='Path folder containing eHive scripts')

