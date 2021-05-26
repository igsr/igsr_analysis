import os
import pytest
import glob
import pdb

from VCF.VcfNormalize import VcfNormalize

# test_VcfNormalize.py

@pytest.fixture
def vcf_object(bcftools_folder, bgzip_folder, gatk_folder, vt_folder, vcflib_folder, datadir):
    """Returns a VcfNormalize object"""

    vcf_file = "{0}/test.vcf.gz".format(datadir)

    vcf_object = VcfNormalize(vcf=vcf_file, bgzip_folder=bgzip_folder,
                              vcflib_folder=vcflib_folder, vt_folder=vt_folder,
                              bcftools_folder=bcftools_folder, gatk_folder=gatk_folder)
    return vcf_object

def test_run_vtnormalize(vcf_object, datadir, clean_tmp):
    outfile = vcf_object.run_vtnormalize(outprefix="{0}/outdir/test".format(datadir),
                                         reference="{0}/exampleFASTA.fasta".format(datadir),
                                         verbose=True,
                                         n=True)

    assert os.path.exists(outfile)

def test_run_vtnormalize_compress(vcf_object, datadir, clean_tmp):
    """
    Test vt normalize with compress=True option
    """
    outfile = vcf_object.run_vtnormalize(outprefix="{0}/outdir/test".format(datadir),
                                         reference="{0}/exampleFASTA.fasta".format(datadir),
                                         verbose=True,
                                         n=True,
                                         compress=True)

    assert os.path.exists(outfile)

def test_run_vcfallelicprimitives(vcf_object, datadir, clean_tmp):

    outfile = vcf_object.run_vcfallelicprimitives(outprefix="{0}/outdir/test".format(datadir),
                                                  verbose=True)

    assert os.path.exists(outfile)

def test_run_vcfallelicprimitives_downstream_pipe(vcf_object, datadir, vt_folder, clean_tmp):
    """
    Test function to run vcfallelicprimitives and piping to other programs
    """
    outfile = vcf_object.run_vcfallelicprimitives(outprefix="{0}/outdir/test1".format(datadir),
                                                  downstream_pipe="{0}/vt sort -".format(vt_folder),
                                                  verbose=True, compress=True)
    assert os.path.exists(outfile)

def test_run_GATK_VariantsToAllelicPrimitives(vcf_object, datadir, clean_tmp):
    """
    Test function to run GATK VariantsToAllelicPrimitives
    """

    outfile = vcf_object.run_gatk_VariantsToAllelicPrimitives(outprefix="{0}/outdir/test2".format(datadir),
                                                              reference="{0}/exampleFASTA.fasta".format(datadir),
                                                              compress=True,
                                                              verbose=True)
    assert os.path.exists(outfile)
