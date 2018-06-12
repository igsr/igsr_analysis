import os
import pytest
import glob
import pdb
import time

from VariantCalling import GATK 

# test_GATK_UG.py

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
    files = glob.glob('data/out/*')
    for f in files:
        os.remove(f)

def test_run_ug(gatk_object):
    '''
    Test function to run GATK UG on a BAM file using a log file
    '''
    
    #create timestamp for log file:
    timestr = time.strftime("%Y%m%d_%H%M%S")

    outfile=gatk_object.run_ug(outprefix='data/out/test', log_file='data/out/gatk_ug_{0}.log'.format(timestr))

    assert os.path.isfile(outfile) is True

def test_run_ug_nocompress(gatk_object):
    '''
    Test function to run GATK UG on a BAM file generating an uncompressed VCF
    '''

    outfile=gatk_object.run_ug(outprefix='data/out/test',verbose=True, compress=False)

    assert os.path.isfile(outfile) is True

def test_run_ug_with_ivals(gatk_object):
    '''
    Test function to run GATK UG on a BAM file and set the interval to 
    analyze on the command line
    '''

    outfile=gatk_object.run_ug(outprefix='data/out/test1', verbose=True, intervals= 'chr1:10000-30000')

    assert os.path.isfile(outfile) is True

def test_run_ug_with_params(gatk_object):
    '''
    Test function to run GATK UG on a BAM file using some optional params
    '''

    outfile=gatk_object.run_ug(outprefix='data/out/test2', glm='INDEL', 
                               output_mode='EMIT_ALL_SITES', nt=1)

    assert os.path.isfile(outfile) is True

def test_run_ug_with_verbose(gatk_object):
    '''
    Test function to run GATK UG on a BAM file using verbose=True
    '''

    outfile=gatk_object.run_ug(outprefix='data/out/test2', glm='INDEL',
                               output_mode='EMIT_ALL_SITES', nt=1, verbose=True)

    assert os.path.isfile(outfile) is True

def test_run_ug_and_throwerror(gatk_object, clean_tmp):
    '''
    Test function to run GATK UG on a BAM file and will raise an Exception 
    because the output_mode argument is not valid
    '''

    #create timestamp for log file:
    timestr = time.strftime("%Y%m%d_%H%M%S")

    with pytest.raises(Exception):
        outfile=gatk_object.run_ug(outprefix='data/out/test2', glm='INDEL',
                                   log_file='data/out/gatk_ug_{0}.log'.format(timestr),
                                   output_mode='non_valid', nt=1)

