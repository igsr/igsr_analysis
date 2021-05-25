import os
import glob
import time
import pytest

from VCF.VCFfilter.GATK import GATK

# test_VcfFilter_GATK.py

@pytest.fixture
def vcf_object(datadir, gatk_folder):
    """Returns a GATK object"""

    vcf_file = "{0}/test.vcf.gz".format(datadir)
    reference = "{0}/exampleFASTA.fasta".format(datadir)

    vcf_object = GATK(vcf=vcf_file, gatk_folder=gatk_folder, reference=reference)
    
    return vcf_object

def test_run_variantrecalibrator(vcf_object, datadir, clean_tmp):
    """
    Test method to run VariantRecalibrator. This test throw an error because the paths
    specified in data/resources_snps.json are fictitious
    """

    #create timestamp for log file:
    timestr = time.strftime("%Y%m%d_%H%M%S")

    with pytest.raises(Exception):
        vcf_object.run_variantrecalibrator(outprefix="{0}/outdir/test".format(datadir),
                                           resources="{0}/resources_snps.json".format(datadir),
                                           mode='SNP', verbose=True,
                                           log_file="{0}/outdir/gatk_variantrecalibrator_{0}.log".format(datadir).
                                           format(timestr))

def test_applyrecalibration(vcf_object, datadir, clean_tmp):
    """
    Test method to run ApplyRecalibration. This test throw an error because the
    'recal_file' and 'tranches_file' files are fictitious
    """

    #create timestamp for log file:
    timestr = time.strftime("%Y%m%d_%H%M%S")

    with pytest.raises(Exception):
        vcf_object.run_applyrecalibration(mode='SNP',
                                          recal_file="{0}/test.recal".format(datadir),
                                          tranches_file="{0}/test.tranches".format(datadir),
                                          outprefix="{0}/outdir/test".format(datadir),
                                          verbose=True,
                                          log_file="{0}/outdir/gatk_applyrecalibration_{0}.log".format(datadir).
                                          format(timestr))

def test_applyrecalibration_uncompressed(vcf_object, datadir, clean_tmp):
    """
    Test method to run ApplyRecalibration in order to generate an uncompressed VCF.
    This test throw an error because the 'recal_file' and 'tranches_file' files are fictitious
    """

    with pytest.raises(Exception):
        vcf_object.run_applyrecalibration(mode='SNP',
                                          recal_file="{0}/test.recal".format(datadir),
                                          tranches_file="{0}/test.tranches".format(datadir),
                                          outprefix="{0}/outdir/test".format(datadir),
                                          compress=False,
                                          verbose=True)

def test_applyrecalibration_tmpdir(datadir, gatk_folder, clean_tmp):
    """
    Test method to run ApplyRecalibration by setting the tmp_dir for Java.
    This test throw an error because the 'recal_file' and 'tranches_file' files are fictitious
    """

    vcf_file = "{0}/test.vcf.gz".format(datadir)
    reference = "{0}/exampleFASTA.fasta".format(datadir)

    vcf_object = GATK(vcf=vcf_file,
                      gatk_folder=gatk_folder,
                      reference=reference,
                      tmp_dir="{0}/tmp".format(datadir))

    with pytest.raises(Exception):
        vcf_object.run_applyrecalibration(mode='SNP',
                                          recal_file="{0}/test.recal".format(datadir),
                                          tranches_file="{0}/test.tranches".format(datadir),
                                          outprefix="{0}/outdir/test".format(datadir),
                                          compress=False,
                                          verbose=True)
