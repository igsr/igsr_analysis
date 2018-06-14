import os
import pytest
import glob
import pdb
import time

from VariantCalling import GATK 

# test_GATK_HC.py

@pytest.fixture
def gatk_object():
    '''Returns an  object'''
    bam_file = pytest.config.getoption("--bam")
    reference = pytest.config.getoption("--reference")
    gatk_folder = pytest.config.getoption("--gatk_folder")
    bgzip_folder = pytest.config.getoption("--bgzip_folder")

    gatk_object=GATK(bam=bam_file, reference=reference, gatk_folder=gatk_folder,bgzip_folder=bgzip_folder)

    return gatk_object

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/outdir/*')
    for f in files:
        os.remove(f)

def test_run_hc(gatk_object):
    '''
    Test function to run GATK HC on a BAM file using a log file
    '''
    
    #create timestamp for log file:
    timestr = time.strftime("%Y%m%d_%H%M%S")

    outfile=gatk_object.run_hc(outprefix='data/outdir/test_hc1', verbose=True, log_file='data/outdir/gatk_hc_{0}.log'.format(timestr))

    assert os.path.isfile(outfile) is True

def test_run_hc_nocompress(gatk_object, clean_tmp):
    '''
    Test function to run GATK HC on a BAM file generating an uncompressed VCF
    '''

    outfile=gatk_object.run_ug(outprefix='data/outdir/test_hc2',verbose=True, compress=False)

    assert os.path.isfile(outfile) is True
