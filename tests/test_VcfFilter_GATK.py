import os
import glob
import time
import pytest

from VCF.VCFfilter.GATK import GATK

# test_VcfFilter_GATK.py

@pytest.fixture
def vcf_object(datadir):
    """Returns a GATK object"""

    vcf_file = "{0}/test.vcf.gz".format(datadir)

    vcf_object = GATK(vcf=vcf_file)
    
    return vcf_object

def test_run_variantrecalibrator(vcf_object, datadir, clean_tmp):
    """
    Test method to run VariantRecalibrator.
    """

    #create timestamp for log file:
    timestr = time.strftime("%Y%m%d_%H%M%S")

    with pytest.raises(Exception):
        vcf_object.run_variantrecalibrator(outprefix="{0}/outdir/test".format(datadir),
                                           resources="{0}/resources_snps.json".format(datadir),
                                           R="{0}/exampleFASTA.fasta".format(datadir),
                                           mode='SNP', verbose=True,
                                           log_file="{0}/outdir/gatk_variantrecalibrator_{0}.log".format(datadir).format(timestr))

def test_applyrecalibration(vcf_object, datadir, clean_tmp):
    """
    Test method to run ApplyRecalibration. This test throw an error because the
    'recal_file' and 'tranches_file' files are fictitious
    """

    #create timestamp for log file:
    timestr = time.strftime("%Y%m%d_%H%M%S")

    with pytest.raises(Exception):
        vcf_object.run_applyrecalibration(outprefix="{0}/outdir/test.filt".format(datadir),
                                          recalFile="{0}/test.recal".format(datadir),
                                          tranchesFile="{0}/test.tranches".format(datadir),
                                          log_file="{0}/outdir/gatk_applyrecalibration_{1}.log".format(datadir,timestr))

def test_applyrecalibration_uncompressed(vcf_object, datadir, clean_tmp):
    """
    Test method to run ApplyRecalibration in order to generate an uncompressed VCF.
    This test throw an error because the 'recal_file' and 'tranches_file' files are fictitious
    """

    #create timestamp for log file:
    timestr = time.strftime("%Y%m%d_%H%M%S")

    with pytest.raises(Exception):
        vcf_object.run_applyrecalibration(outprefix="{0}/outdir/test.filt".format(datadir),
                                          recalFile="{0}/test.recal".format(datadir),
                                          tranchesFile="{0}/test.tranches".format(datadir),
                                          compress=False,
                                          log_file="{0}/outdir/gatk_applyrecalibration_{1}.log".format(datadir,timestr))
