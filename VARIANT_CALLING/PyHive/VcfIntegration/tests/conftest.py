import pytest

def pytest_addoption(parser):
    parser.addoption('--vcf', default='data/test_chr20.vcf.gz' ,action='store_true', help='Path to vcf file')
    parser.addoption('--beagle_folder', default='~/bin/beagle/' ,action='store_true', help='Folder with Beagle jar file')
    parser.addoption('--prepareGenFromBeagle4_folder', default='~/bin/shapeit2_v2_12/bin/prepareGenFromBeagle4/bin/', action='store_true', help='Folder with prepareGenFromBeagle4 binary')
    parser.addoption('--hive_lib', default='~/lib/ensembl-hive_2.4/' ,action='store_true', help='Path folder containing eHive scripts')



