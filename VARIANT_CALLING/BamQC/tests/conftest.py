import pytest

def pytest_addoption(parser):
    parser.addoption('--bam', default='data/exampleBAM.bam' ,action='store_true', help='Path to test BAM file')
    parser.addoption('--chk_indel_folder', default='~/bin/' ,action='store_true', help='Folder with chk_indel_rg binary')
    parser.addoption('--samtools_folder', default='/Users/ernesto/bin/anaconda/envs/python3/bin/',
                     action='store_true', help='Folder with samtools binary')
    parser.addoption('--java_folder', default='/usr/bin/',action='store_true', help='Folder with java binary')
    parser.addoption('--picard_folder', default='~/bin/',action='store_true', help='Folder with Picard jar file')



