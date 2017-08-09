import os
import pytest
from VcfUtils import VcfUtils

# test_VcfUtils.py

@pytest.fixture
def vcf_object():
    '''Returns an  object'''
    vcf_file = pytest.config.getoption("--vcf")
    bcftools_folder = pytest.config.getoption("--bcftools_folder")
    vcf_object=VcfUtils(vcf=vcf_file,bcftools_folder=bcftools_folder)
    return vcf_object


@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    os.remove('data/test.reheaded.vcf.gz')


def test_vcf_reheader(vcf_object, clean_tmp):
    outfile=vcf_object.reheader(newheader="data/newheader.txt", outprefix="data/test") 

    assert os.path.exists("data/test.reheaded.vcf.gz")

def test_vcf_reheader_with_samplef(vcf_object, clean_tmp):
    '''
    Test the reheader method and add new sample names
    '''
    outfile=vcf_object.reheader(newheader="data/newheader.txt", samplefile="data/samples.txt",
                                outprefix="data/test")

    assert os.path.exists("data/test.reheaded.vcf.gz")
