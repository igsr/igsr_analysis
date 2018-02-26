import os
import pytest
import glob
import pdb

from GATK import GATK 

# test_GATK.py

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
    files = glob.glob('out/*')
    for f in files:
        os.remove(f)

def test_run_ug(gatk_object):
    '''
    Test function to run GATK UG on a BAM file
    '''
    
    outfile=gatk_object.run_ug(outprefix='./out/test')

    assert os.path.isfile(outfile) is True

def test_run_ug_with_ivals(gatk_object):
    '''
    Test function to run GATK UG on a BAM file and set the interval to 
    analyze on the command line
    '''

    outfile=gatk_object.run_ug(outprefix='./out/test1', intervals= 'chr1:10000-30000')

    assert os.path.isfile(outfile) is True

def test_run_ug_with_params(gatk_object, clean_tmp):
    '''
    Test function to run GATK UG on a BAM file using some optional params
    '''

    outfile=gatk_object.run_ug(outprefix='./out/test2', glm='INDEL', 
                               output_mode='EMIT_ALL_SITES', num_threads=1)

    assert os.path.isfile(outfile) is True

