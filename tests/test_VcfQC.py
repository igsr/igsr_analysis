import os
import pytest
from VcfQC import VcfQC

# test_vcfqc.py

@pytest.fixture
def vcf_object():
    '''Returns an  object'''
    vcf_file = pytest.config.getoption("--vcf")
    bcftools_folder = pytest.config.getoption("--bcftools_folder")
    picard_folder = pytest.config.getoption("--picard_folder")
    vcf_object=VcfQC(vcf=vcf_file,bcftools_folder=bcftools_folder,picard_folder=picard_folder)
    return vcf_object

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    os.remove('data/outdir/test_VcfQC.stats')

def test_stats(vcf_object):
    '''
    Test function to run Bcftools stats on a VCF file
    '''
    stats=vcf_object.stats(outpath='data/outdir/test_VcfQC', verbose=True) 
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

def test_stats_with_region(vcf_object):
    '''
    Test function to run Bcftools stats on a VCF file
    for a certain genomic region
    '''
    stats=vcf_object.stats(outpath='data/outdir/test_VcfQC', region='chr1:80000-90000', verbose=True)
    assert os.path.isfile(stats.filename) is True
    assert stats.summary_numbers['number of SNPs:']==9
    assert stats.summary_numbers['number of indels:']==3
    
def test_get_chros(vcf_object):
    chr_file = pytest.config.getoption("--chr_file")
    dict=vcf_object.get_chros(chr_f=chr_file, verbose=True)
    assert dict['both'][0] == "chr1"

def test_number_variants_in_region(vcf_object):
    counts=vcf_object.number_variants_in_region(outprefix='data/outdir/test_number_vars_in_region',region= 
                                               pytest.config.getoption("--region"), verbose=True)
    
    assert os.path.exists(counts)

def test_run_CollectVariantCallingMetrics(vcf_object,clean_tmp):
    '''
    Test function for running Picard CollectVariantCallingMetrics
    '''
    
    cvcm=vcf_object.run_CollectVariantCallingMetrics(outprefix='data/outdir/test', truth_vcf=vcf_object.vcf, verbose=True)

    assert cvcm.vc_detail_metrics['DBSNP_TITV']=='2.666667'
    assert cvcm.vc_detail_metrics['HET_HOMVAR_RATIO']=='2.84375'


