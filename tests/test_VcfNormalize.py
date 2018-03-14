import os
import pytest
import glob
import pdb

from VcfNormalize import VcfNormalize

# test_VcfNormalize.py

@pytest.fixture
def vcf_object():
    '''Returns an  object'''

    vcf_file = pytest.config.getoption("--vcf")
    bcftools_folder = pytest.config.getoption("--bcftools_folder")
    bgzip_folder = pytest.config.getoption("--bgzip_folder")
    vcflib_folder = pytest.config.getoption("--vcflib_folder")

    vcf_object=VcfNormalize(vcf=vcf_file,bgzip_folder=bgzip_folder,
                            vcflib_folder=vcflib_folder, bcftools_folder=bcftools_folder)
    return vcf_object

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('out/*')
    for f in files:
        os.remove(f)

def test_run_vcfallelicprimitives(vcf_object):
    outfile=vcf_object.run_vcfallelicprimitives(outprefix='out/test')

    assert os.path.exists("out/test.aprimitives.vcf")

def test_run_vcfallelicprimitives_downstream_pipe(vcf_object):
    '''
    Test funtion to run vcfallelicprimitives and piping to other programs
    '''
    outfile=vcf_object.run_vcfallelicprimitives(outprefix='out/test1', 
                                                downstream_pipe='~/bin/vt/vt sort -',
                                                compress=True)

    assert os.path.exists("out/test1.aprimitives.vcf")
