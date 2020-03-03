import os
import glob
import pytest
import pdb

from VariantCalling import BCFTools

# test_BCFTools.py

@pytest.fixture
def bcftools_object(supply_bam_file, supply_reference_file, supply_settings_file):
    '''Returns a BCFTools object'''
    bcftools_object = BCFTools(bam=supply_bam_file,
                               reference=supply_reference_file,
                               settings=supply_settings_file)

    return bcftools_object

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/outdir/*')
    for f in files:
        os.remove(f)

def test_run_bcftools(bcftools_object, clean_tmp):
    '''
    Test function to run BCFTools on a BAM file in order to call variants
    '''

    outfile = bcftools_object.run_bcftools(outprefix='data/outdir/test')

    assert os.path.isfile(outfile) is True
"""
def test_run_bcftools_w_region(bcftools_object, clean_tmp):
    '''
    Test function to run BCFTools on a BAM file for a certain region
    '''

    annots = ['DP', 'SP', 'AD']

    outfile = bcftools_object.run_bcftools(outprefix='data/outdir/test',
                                           annots=annots,
                                           E=True,
                                           p=True,
                                           m_pileup=3,
                                           m_call=True,
                                           r="chr1:400-1000",
                                           v=True)

    assert os.path.isfile(outfile) is True

def test_run_bcftools_w_threads(bcftools_object, clean_tmp):
    '''
    Test function to run BCFTools on a BAM file using more than one thread
    '''

    annots = ['DP', 'SP', 'AD']

    outfile = bcftools_object.run_bcftools(outprefix='data/outdir/test',
                                           annots=annots,
                                           E=True,
                                           p=True,
                                           m_pileup=3,
                                           m_call=True,
                                           r="chr1:400-1000",
                                           threads=2,
                                           v=True)

    assert os.path.isfile(outfile) is True
"""
