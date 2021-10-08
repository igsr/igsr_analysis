import os
import time
import pytest
import pdb

from VariantCalling import GATK

# test_GATK_HC.py

def test_run_hc(gatk_object, datadir, clean_tmp):
    """
    Test function to run GATK HC on a BAM file using a log file
    """
    #create timestamp for log file:
    timestr = time.strftime("%Y%m%d_%H%M%S")

    outfile = gatk_object.run_caller(program='HaplotypeCaller',
                                     prefix="{0}/outdir/test_hc1".format(datadir),
                                     log_file="{0}/outdir/gatk_hc_{1}.log".format(datadir, timestr))

    assert os.path.isfile(outfile) is True

def test_run_hc_nocompress(gatk_object, datadir, clean_tmp):
    """
    Test function to run GATK HC on a BAM file generating an uncompressed VCF
    """

    outfile = gatk_object.run_caller(program='HaplotypeCaller',
                                     prefix="{0}/outdir/test_hc2".format(datadir),
                                     compress=False)

    assert os.path.isfile(outfile) is True

def test_run_hc_interval(gatk_object, datadir, clean_tmp):
    """
    Test function to run GATK HC on a BAM file using several intervals
    """

    outfile = gatk_object.run_caller(program='HaplotypeCaller',
                                     prefix="{0}/outdir/test_hc3".format(datadir),
                                     L='chr1:10000-20000')

    assert os.path.isfile(outfile) is True
