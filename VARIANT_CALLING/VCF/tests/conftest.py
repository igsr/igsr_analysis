import pytest

def pytest_addoption(parser):
    parser.addoption('--vcf', default='data/test.vcf.gz' ,action='store_true', help='Path to test VCF file')
    parser.addoption('--bcftools_folder', default='~/bin/bcftools/' ,action='store_true', help='Folder with bcftools binary')
    parser.addoption('--snptools_folder', default='~/bin/snptools/' ,action='store_true', help='Folder with SNPTools binaries')
    parser.addoption('--beagle_folder', default='~/bin/beagle/' ,action='store_true', help='Folder with Beagle jar file')
    parser.addoption('--chr_file', default='data/chr_file.txt' ,action='store_true', help='File with chros for get_chros function')

