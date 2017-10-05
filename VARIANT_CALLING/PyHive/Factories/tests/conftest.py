import pytest

def pytest_addoption(parser):
    parser.addoption('--faix', default='data/canonical_chros.fa.fai' ,action='store_true', help='Path to faix file')
    parser.addoption('--vcf_gts', default='data/GLs.HG00136.vcf.gz' ,action='store_true', help='Path to vcf file with GTs')
    parser.addoption('--vcf_gts_ucsc', default='data/GLs.HG00136.ucsc.vcf.gz' ,action='store_true', help='Path to vcf file with GTs with UCSC-style chro names')
    parser.addoption('--makeBGLCHUNKS_folder', default='~/bin/shapeit2_v2_12/bin/makeBGLCHUNKS/bin/' ,action='store_true', help='Folder with makeBGLCHUNKS binary')
    parser.addoption('--bedtools_folder', default='/homes/ernesto/bin/bedtools-2.25.0/bin/', action='store_true', help='Folder with bedtools binary')
    parser.addoption('--hive_lib', default='~/lib/ensembl-hive_2.4/' ,action='store_true', help='Path folder containing eHive scripts')

