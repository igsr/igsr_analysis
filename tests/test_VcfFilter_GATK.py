import os
import pytest
import glob

from VCFfilter import GATK

# test_VcfFilter_GATK.py

@pytest.fixture
def vcf_object():
    '''Returns an  object'''
    vcf_file = pytest.config.getoption("--vcf")
    gatk_folder = pytest.config.getoption("--gatk_folder")
    reference = pytest.config.getoption("--reference")
    vcf_object=GATK(vcf=vcf_file,gatk_folder=gatk_folder,reference=reference)
    return vcf_object

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/out/*')
    for f in files:
        os.remove(f)

def test_run_variantrecalibrator(vcf_object):
    '''
    Test method to run VariantRecalibrator. This test throw an error because the paths specified in data/resources_snps.json 
    are fictitious
    '''

    with pytest.raises(Exception):
        vcf_object.run_variantrecalibrator(outprefix='data/outdir/test', resources='data/resources_snps.json', mode='SNP', verbose=True)

def test_applyrecalibration(vcf_object):
    '''
    Test method to run ApplyRecalibration. This test throw an error because the 'recal_file' and 'tranches_file' files  are ficticious
    '''
    
    with pytest.raises(Exception):
        vcf_object.run_applyrecalibration(mode='SNP', recal_file='data/test.recal', 
                                          tranches_file='data/test.tranches', outprefix='data/outdir/test',
                                          verbose=True)

def test_applyrecalibration_uncompressed(vcf_object):
    '''
    Test method to run ApplyRecalibration in order to generate an uncompressed VCF.
    This test throw an error because the 'recal_file' and 'tranches_file' files  are ficticious
    '''

    with pytest.raises(Exception):
        vcf_object.run_applyrecalibration(mode='SNP', recal_file='data/test.recal',
                                          tranches_file='data/test.tranches', outprefix='data/outdir/test',
                                          compress=False, verbose=True)    

