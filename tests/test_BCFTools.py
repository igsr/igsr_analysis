import os
import glob
import pytest
import shutil
import pdb

from VariantCalling import BCFTools

# test_BCFTools.py

@pytest.fixture
def bcftools_object(datadir, bcftools_folder):
    '''Returns a BCFTools object'''

    bam_file = "{0}/exampleBAM.bam".format(datadir)
    reference = "{0}/exampleFASTA.fasta".format(datadir)

    bcftools_object = BCFTools(bam=bam_file, reference=reference, bcftools_folder=bcftools_folder)

    return bcftools_object

@pytest.fixture
def clean_tmp(datadir):
    yield
    print("Cleanup files")
    files = glob.glob("{0}/outdir/*".format(datadir))
    for f in files:
        os.remove(f)

def test_run_bcftools(bcftools_object, datadir, clean_tmp):
    '''
    Test function to run BCFTools on a BAM file in order to call variants
    '''

    annots = ['DP', 'SP', 'AD']

    outfile = bcftools_object.run_bcftools(outprefix="{0}/outdir/test".format(datadir),
                                           annots=annots,
                                           E=True,
                                           p=True,
                                           m_pileup=3,
                                           m_call=True,
                                           v=True)

    assert os.path.isfile(outfile) is True

def test_run_bcftools_w_region(bcftools_object, datadir, clean_tmp):
    '''
    Test function to run BCFTools on a BAM file for a certain region
    '''

    annots = ['DP', 'SP', 'AD']

    outfile = bcftools_object.run_bcftools(outprefix="{0}/outdir/test".format(datadir),
                                           annots=annots,
                                           E=True,
                                           p=True,
                                           m_pileup=3,
                                           m_call=True,
                                           r="chr1:400-1000",
                                           v=True)

    assert os.path.isfile(outfile) is True

def test_run_bcftools_w_threads(bcftools_object, datadir, clean_tmp):
    '''
    Test function to run BCFTools on a BAM file using more than one thread
    '''

    annots = ['DP', 'SP', 'AD']

    outfile = bcftools_object.run_bcftools("{0}/outdir/test".format(datadir),
                                           annots=annots,
                                           E=True,
                                           p=True,
                                           m_pileup=3,
                                           m_call=True,
                                           r="chr1:400-1000",
                                           threads=2,
                                           v=True)

    assert os.path.isfile(outfile) is True
