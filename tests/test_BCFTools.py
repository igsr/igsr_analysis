import os
import glob
import pytest
import shutil
import pdb

from VariantCalling import BCFTools

# test_BCFTools.py

@pytest.fixture
def bcftools_object(datadir):
    '''Returns a BCFTools object'''

    bam_file = "{0}/exampleBAM.bam".format(datadir)
    reference = "{0}/exampleFASTA.fasta".format(datadir)

    bcftools_object = BCFTools(bam=bam_file, reference=reference)

    return bcftools_object

def test_run_bcftools(bcftools_object, datadir, clean_tmp):
    '''
    Test function to run BCFTools on a BAM file in order to call variants
    '''

    annots = ['DP', 'SP', 'AD']

    outfile = bcftools_object.run_bcftools('m',
                                           'E',
                                           prefix="{0}/outdir/test".format(datadir),
                                           annots=annots)

    assert os.path.isfile(outfile) is True

def test_run_bcftools_w_region(bcftools_object, datadir, clean_tmp):
    '''
    Test function to run BCFTools on a BAM file for a certain region
    '''

    annots = ['DP', 'SP', 'AD']

    outfile = bcftools_object.run_bcftools('m',
                                           'e',
                                           prefix="{0}/outdir/test".format(datadir),
                                           annots=annots,
                                           r="chr1:400-1000")

    assert os.path.isfile(outfile) is True

def test_run_bcftools_w_threads(bcftools_object, datadir, clean_tmp):
    '''
    Test function to run BCFTools on a BAM file using more than one thread
    '''

    annots = ['DP', 'SP', 'AD']

    outfile = bcftools_object.run_bcftools('m',
                                           'e',
                                           prefix="{0}/outdir/test".format(datadir),
                                           annots=annots,
                                           r="chr1:400-1000",
                                           threads=2)

    assert os.path.isfile(outfile) is True
