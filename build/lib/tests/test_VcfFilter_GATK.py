import os
import pytest
import glob
import time

from VCF.VCFfilter.GATK import GATK

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
    files = glob.glob('data/outdir/*')
    for f in files:
        os.remove(f)

def test_run_variantrecalibrator(vcf_object):
    '''
    Test method to run VariantRecalibrator. This test throw an error because the paths specified in data/resources_snps.json 
    are fictitious
    '''

    #create timestamp for log file:
    timestr = time.strftime("%Y%m%d_%H%M%S")

    with pytest.raises(Exception):
        vcf_object.run_variantrecalibrator(outprefix='data/outdir/test', resources='data/resources_snps.json', mode='SNP', verbose=True,
                                           log_file='data/outdir/gatk_variantrecalibrator_{0}.log'.format(timestr))

def test_applyrecalibration(vcf_object):
    '''
    Test method to run ApplyRecalibration. This test throw an error because the 'recal_file' and 'tranches_file' files  are ficticious
    '''

    #create timestamp for log file:
    timestr = time.strftime("%Y%m%d_%H%M%S")
    
    with pytest.raises(Exception):
        vcf_object.run_applyrecalibration(mode='SNP', recal_file='data/test.recal', 
                                          tranches_file='data/test.tranches', outprefix='data/outdir/test',
                                          verbose=True,log_file='data/outdir/gatk_applyrecalibration_{0}.log'.format(timestr))

def test_applyrecalibration_uncompressed(vcf_object):
    '''
    Test method to run ApplyRecalibration in order to generate an uncompressed VCF.
    This test throw an error because the 'recal_file' and 'tranches_file' files  are ficticious
    '''

    with pytest.raises(Exception):
        vcf_object.run_applyrecalibration(mode='SNP', recal_file='data/test.recal',
                                          tranches_file='data/test.tranches', outprefix='data/outdir/test',
                                          compress=False, verbose=True)

def test_applyrecalibration_tmpdir(clean_tmp):
    '''
    Test method to run ApplyRecalibration by setting the tmp_dir for Java.
    This test throw an error because the 'recal_file' and 'tranches_file' files  are ficticious
    '''

    vcf_file = pytest.config.getoption("--vcf")
    gatk_folder = pytest.config.getoption("--gatk_folder")
    reference = pytest.config.getoption("--reference")
    vcf_object=GATK(vcf=vcf_file,gatk_folder=gatk_folder,reference=reference,tmp_dir='data/tmp')

    with pytest.raises(Exception):
        vcf_object.run_applyrecalibration(mode='SNP', recal_file='data/test.recal',
                                          tranches_file='data/test.tranches', outprefix='data/outdir/test',
                                          compress=False, verbose=True)


