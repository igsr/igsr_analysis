import os
import glob
import pytest

from VCF.VCFfilter.BCFTools import BCFTools

# test_vcfFilter_BCFTools.py

@pytest.fixture
def vcf_object(datadir):
    """Returns a BCFTools object"""

    vcf_file = "{0}/test.vcf.gz".format(datadir)
    vcf_object = BCFTools(vcf=vcf_file)

    return vcf_object

def test_filter_by_variant_type(vcf_object, datadir, clean_tmp):
    """
    Test method filter_by_variant_type
    Will select SNPs from the VCF file
    """

    outfile = vcf_object.filter_by_variant_type(prefix="{0}/outdir/test".format(datadir))

    assert os.path.isfile(outfile) is True

def test_filter_by_variant_type_biallelic(vcf_object, datadir, clean_tmp):
    """
    Test method filter_by_variant_type
    using the biallelic option
    """

    outfile = vcf_object.filter_by_variant_type(prefix="{0}/outdir/test".format(datadir),
                                                biallelic=True,
                                                verbose=True)

    assert os.path.isfile(outfile) is True

def test_filter_by_variant_type_biallelic_compressed(vcf_object, datadir, clean_tmp):
    """
    Test method filter_by_variant_type
    using the biallelic option
    """

    outfile = vcf_object.filter_by_variant_type(prefix="{0}/outdir/test".format(datadir),
                                                biallelic=True,
                                                compress=False,
                                                verbose=True)

    assert os.path.isfile(outfile) is True

def test_subset_vcf(vcf_object, datadir, clean_tmp):
    """
    Test method subset_vcf to subset a VCF by using a BED file/region
    """

    outfile = vcf_object.subset_vcf(prefix="{0}/outdir/test.vcf.gz".format(datadir),
                                    r="chr1")

    assert os.path.isfile(outfile) is True

def test_subset_vcf_and_throwerror(vcf_object, datadir, clean_tmp):
    """
    Test method subset_vcf to subset a VCF by using a BED file/region and using an invalid 'action'
    parameter to throw an exception
    """

    with pytest.raises(Exception):
        outfile = vcf_object.subset_vcf(outprefix="{0}/outdir/test.vcf.gz".format(datadir),
                                        region="chr1",
                                        action='test',
                                        apply_filters="PASS",
                                        verbose=True)

def test_select_variants(vcf_object, datadir, clean_tmp):
    """
    Test method to select only the variants (exclude the 0|0 genotypes) from a VCF file
    """

    outfile = vcf_object.select_variants(outprefix="{0}/outdir/test".format(datadir))

    assert os.path.isfile(outfile) is True

def test_select_variants_exclude_uncalled(vcf_object, datadir, clean_tmp):
    """
    Test method to select only the variants (exclude the 0|0 genotypes) and also exclude
    the sites with uncalled genoytpes from a VCF file
    """

    outfile = vcf_object.select_variants(outprefix="{0}/outdir/test".format(datadir),
                                         uncalled='exclude',
                                         verbose=True)

    assert os.path.isfile(outfile) is True

def test_filter(vcf_object, clean_tmp):
    """
    Test method to filter variants from a VCF file by running bcftools filter
    """

    outfile = vcf_object.filter(name='TESTFILTER', expression="'INFO/DP>24304'")

    assert os.path.isfile(outfile) is True
