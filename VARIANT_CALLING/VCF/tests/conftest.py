import pytest

def pytest_addoption(parser):
    parser.addoption('--vcf', default='data/test.vcf.gz', action='store_true', help='Path to test VCF file')
    parser.addoption('--vcf_ambiguity', default='data/test.amb.vcf.gz', action='store_true', help='Path to VCF file containing the REF or ALT column with some ambiguity codes')
    parser.addoption('--vcflist', default=['data/test.vcf.gz','data/test1.vcf.gz'], action='store_true', help='List with VCF paths')
    parser.addoption('--bcftools_folder', default='/nfs/software/ensembl/RHEL7-JUL2017-core2/linuxbrew/bin/', action='store_true', help='Folder with bcftools binary')
    parser.addoption('--gatk_folder', default='~/bin/GATK/', action='store_true', help='Folder with GATK jar file')
    parser.addoption('--snptools_folder', default='~/bin/snptools/', action='store_true', help='Folder with SNPTools binaries')
    parser.addoption('--beagle_folder', default='~/bin/beagle/', action='store_true', help='Folder with Beagle jar file')
    parser.addoption('--makeBGLCHUNKS_folder', default='~/bin/shapeit2_v2_12/bin/makeBGLCHUNKS/bin/', action='store_true', help='Folder with makeBGLCHUNKS binary')
    parser.addoption('--prepareGenFromBeagle4_folder', default='~/bin/shapeit2_v2_12/bin/prepareGenFromBeagle4/bin/', action='store_true', help='Folder with prepareGenFromBeagle4 binary')
    parser.addoption('--shapeit_folder', default='~/bin/shapeit2_v2_12/bin/', action='store_true', help='Folder with Shapeit binary')
    parser.addoption('--chr_file', default='data/chr_file.txt', action='store_true', help='File with chros for get_chros function')

