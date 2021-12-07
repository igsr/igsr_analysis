import os
import time
import pytest
import pdb

from VariantCalling import GATK

# test_GATK_UG.py

def test_run_ug(gatk_object, datadir, clean_tmp):
    """
    Test function to run GATK UG on a BAM file using a log file
    """

    #create timestamp for log file:
    timestr = time.strftime("%Y%m%d_%H%M%S")

    outfile = gatk_object.run_caller(program='UnifiedGenotyper',
                                     prefix="{0}/outdir/test".format(datadir),
                                     log_file="{0}/outdir/gatk_ug_{1}.log".format(datadir, timestr))

    assert os.path.isfile(outfile) is True

def test_run_ug_nocompress(gatk_object, datadir, clean_tmp):
    """
    Test function to run GATK UG on a BAM file generating an uncompressed VCF
    """

    outfile = gatk_object.run_caller(program='UnifiedGenotyper', 
                                 prefix="{0}/outdir/test".format(datadir),
                                 compress=False)

    assert os.path.isfile(outfile) is True

def test_run_ug_with_ivals(gatk_object, datadir, clean_tmp):
    """
    Test function to run GATK UG on a BAM file and set the interval to
    analyze on the command line
    """

    outfile = gatk_object.run_caller(program='UnifiedGenotyper',
                                     prefix="{0}/outdir/test1".format(datadir),
                                     L='chr1:10000-20000')

    assert os.path.isfile(outfile) is True

def test_run_ug_with_params(gatk_object, datadir, clean_tmp):
    """
    Test function to run GATK UG on a BAM file using some optional params
    """

    outfile = gatk_object.run_caller(program='UnifiedGenotyper',
                                     prefix="{0}/outdir/test2".format(datadir), glm='INDEL',
                                     out_mode='EMIT_ALL_SITES')

    assert os.path.isfile(outfile) is True

def test_run_ug_with_verbose(gatk_object, datadir, clean_tmp):
    """
    Test function to run GATK UG on a BAM file using verbose=True
    """

    outfile = gatk_object.run_caller(program='UnifiedGenotyper',
                                     prefix="{0}/outdir/test2".format(datadir),
                                     verbose=True)

    assert os.path.isfile(outfile) is True

def test_run_ug_multithreaded(gatk_object, datadir, clean_tmp):
    """
    Test function to run GATK UG on a BAM file using more than one thread
    """

    outfile = gatk_object.run_caller(program='UnifiedGenotyper',
                                     prefix="{0}/outdir/test2".format(datadir),
                                     nt=2)
