import os
import pytest
import pdb

from VCF.VcfQC import VcfQC

# test_vcfqc.py

@pytest.fixture
def vcf_object(picard_folder, datadir, bcftools_folder):
    """Returns a VcfQC object"""
    
    vcf_file = "{0}/test.vcf.gz".format(datadir)

    vcf_object = VcfQC(vcf=vcf_file, bcftools_folder=bcftools_folder, picard_folder=picard_folder)
    
    return vcf_object

def test_stats(vcf_object, datadir, clean_tmp):
    """
    Test function to run Bcftools stats on a VCF file
    """

    stats = vcf_object.stats(outpath="{0}/outdir/test_VcfQC".format(datadir), verbose=True)
    assert os.path.isfile(stats.filename) is True
    assert stats.summary_numbers['number of SNPs:'] == 112
    assert stats.summary_numbers['number of indels:'] == 13
    assert stats.summary_numbers['number of records:'] == 123
    assert stats.summary_numbers['number of samples:'] == 1
    assert stats.summary_numbers['number of others:'] == 3
    assert stats.summary_numbers['number of no-ALTs:'] == 0
    assert stats.summary_numbers['number of multiallelic SNP sites:'] == 4
    assert stats.summary_numbers['number of MNPs:'] == 10
    assert stats.summary_numbers['number of multiallelic sites:'] == 20

def test_stats_with_region(vcf_object, datadir, clean_tmp):
    """
    Test function to run Bcftools stats on a VCF file
    for a certain genomic region
    """
    stats = vcf_object.stats(outpath="{0}/outdir/test_VcfQC".format(datadir),
                             region='chr1:80000-90000',
                             verbose=True)

    assert os.path.isfile(stats.filename) is True
    assert stats.summary_numbers['number of SNPs:'] == 9
    assert stats.summary_numbers['number of indels:'] == 3

def test_get_chros(vcf_object, datadir, clean_tmp):

    chr_file = "{0}/chr_file.txt".format(datadir)
    dict = vcf_object.get_chros(chr_f=chr_file, verbose=True)
    assert dict['both'][0] == "chr1"

def test_number_variants_in_region(vcf_object, datadir, clean_tmp):

    counts = vcf_object.number_variants_in_region(outprefix="{0}/outdir/test_number_vars_in_region".format(datadir),
                                                  region="{0}/region.bed".format(datadir),
                                                  verbose=True)

    assert os.path.exists(counts)

def test_run_CollectVariantCallingMetrics(vcf_object, datadir, clean_tmp):
    """
    Test function for running Picard CollectVariantCallingMetrics
    """

    cvcm = vcf_object.run_CollectVariantCallingMetrics(outprefix="{0}/outdir/test".format(datadir),
                                                       truth_vcf=vcf_object.vcf,
                                                       verbose=True)

    assert cvcm.vc_detail_metrics['DBSNP_TITV'] == '2.666667'
    assert cvcm.vc_detail_metrics['HET_HOMVAR_RATIO'] == '2.84375'
