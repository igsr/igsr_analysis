import os
import pytest
from VcfQC import VcfQC

# test_vcfqc.py

@pytest.fixture
def vcf_object():
    '''Returns an  object'''
    vcf_file = pytest.config.getoption("--vcf")
    bcftools_folder = pytest.config.getoption("--bcftools_folder")
    vcf_object=VcfQC(vcf=vcf_file,bcftools_folder=bcftools_folder)
    return vcf_object

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    os.remove('data/out/test_VcfQC.stats')

def test_bcftools_stats(vcf_object, clean_tmp):
    stats=vcf_object.stats(outpath='data/out/test_VcfQC') 
    assert os.path.isfile(stats.filename) is True
    assert stats.summary_numbers['number of SNPs:']==112
    assert stats.summary_numbers['number of indels:']==13
    assert stats.summary_numbers['number of records:']==123
    assert stats.summary_numbers['number of samples:']==1
    assert stats.summary_numbers['number of others:']==3
    assert stats.summary_numbers['number of no-ALTs:']==0
    assert stats.summary_numbers['number of multiallelic SNP sites:']==4
    assert stats.summary_numbers['number of MNPs:']==10
    assert stats.summary_numbers['number of multiallelic sites:']==20

def test_get_chros(vcf_object):
    chr_file = pytest.config.getoption("--chr_file")
    dict=vcf_object.get_chros(chr_f=chr_file)
    assert dict['both'][0] == "chr1"


