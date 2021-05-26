import os
import glob
import pytest
import pdb

from VCF.VcfUtils import VcfUtils

# test_VcfUtils.py

@pytest.fixture
def vcf_object(bcftools_folder, bgzip_folder, gatk_jar_folder, datadir):
    """Returns a VcfUtils object"""

    vcf_file = "{0}/test.vcf.gz".format(datadir)
    vcflist = ['test.vcf.gz','test1.vcf.gz']
    vcflist = [ "{0}/{1}".format(datadir, x) for x in vcflist]
   
    vcf_object = VcfUtils(vcf=vcf_file,
                          bgzip_folder=bgzip_folder,
                          vcflist=vcflist,
                          bcftools_folder=bcftools_folder,
                          gatk_folder=gatk_jar_folder)
    return vcf_object

@pytest.fixture
def vcf_ambiguity(datadir):
    """
    Returns a VcfUtils object that contains the REF or ALT column with some ambiguity codes
    """
    vcf_file = "{0}/test.amb.vcf.gz".format(datadir)
    vcf_object = VcfUtils(vcf=vcf_file)

    return vcf_object

def test_correct_ambiguity(vcf_ambiguity, datadir):

    outfile = vcf_ambiguity.correct_ambiguity_codes(outfile="{0}/outdir/test.corrected.vcf.gz".format(datadir))

    assert os.path.exists("{0}/outdir/test.corrected.vcf.gz".format(datadir))

def test_add_to_header(vcf_ambiguity, datadir):

    outfile = vcf_ambiguity.add_to_header(header_f="{0}/newheader.txt".format(datadir),
                                          outfilename="{0}/outdir/modified_header.txt".format(datadir),
                                          line_ann='##INFO=test"')

    assert os.path.exists("{0}/outdir/modified_header.txt".format(datadir))

def test_vcf_reheader(vcf_object, datadir):
    outfile = vcf_object.reheader(newheader="{0}/newheader.txt".format(datadir),
                                  outprefix="{0}/outdir/test1".format(datadir),
                                  verbose=True)

    assert os.path.exists("{0}/outdir/test1.reheaded.vcf.gz".format(datadir))

def test_vcf_reheader_with_samplef(vcf_object, datadir):
    """
    Test the reheader method and add new sample names
    """
    outfile = vcf_object.reheader(newheader="{0}/newheader.txt".format(datadir),
                                  samplefile="{0}/samples.txt".format(datadir),
                                  outprefix="{0}/outdir/test2".format(datadir))

    assert os.path.exists("{0}/outdir/test2.reheaded.vcf.gz".format(datadir))

def test_combine_uncompressed(vcf_object, datadir, clean_tmp):
    """
    Test the combine method producing a VCF
    """
    vcf_object.combine(labels=['gatk', 'lc_bcftools'],
                       reference="{0}/exampleFASTA.fasta".format(datadir),
                       outprefix='out_combine',
                       outdir="{0}/outdir/".format(datadir),
                       verbose=True,
                       genotypemergeoption='UNIQUIFY')

    assert os.path.exists("{0}/outdir/out_combine.vcf".format(datadir))
    assert os.path.exists("{0}/outdir/out_combine.vcf.idx".format(datadir))

def test_combine_compressed(vcf_object, datadir, clean_tmp):
    """
    Test the combine method producing a VCF.gz file and passing also
    some options
    """
    vcf_object.combine(labels=['gatk', 'lc_bcftools'],
                       reference="{0}/exampleFASTA.fasta".format(datadir),
                       outprefix='out_combine',
                       outdir="{0}/outdir/".format(datadir),
                       compress=True,
                       genotypemergeoption='UNIQUIFY',
                       options=['-env', '-sites_only', '--filteredAreUncalled'])

    assert os.path.exists("{0}/outdir/out_combine.vcf.gz".format(datadir))

def test_change_chrnames_2ensembl(vcf_object, datadir):
    """
    Test the method to change the style of the chrnames (from UCSC to Ensembl)
    """

    vcf_object.rename_chros(chr_types='ensembl', outfile="{0}/outdir/test.ensembl.vcf.gz".format(datadir))
    vcf_object.rename_chros(chr_types='ensembl',
                            compress=False,
                            outfile="{0}/outdir/test.ensembl.vcf".format(datadir))

    assert os.path.exists("{0}/outdir/test.ensembl.vcf.gz".format(datadir))
    assert os.path.exists("{0}/outdir/test.ensembl.vcf".format(datadir))

def test_change_chrnames_2ucsc(datadir, bgzip_folder, clean_tmp):
    """
    Test the method to change the style of the chrnames (from Ensembl to UCSC)
    """
    vcf_object = VcfUtils(vcf="{0}/outdir/test.ensembl.vcf.gz".format(datadir),
                          bgzip_folder=bgzip_folder)

    vcf_object.rename_chros(chr_types='ucsc', outfile="{0}/outdir/test.ucsc.vcf.gz".format(datadir))
    vcf_object.rename_chros(chr_types='ucsc', outfile="{0}/outdir/test.ucsc.vcf".format(datadir), compress=False)

    assert os.path.exists("{0}/outdir/test.ucsc.vcf.gz".format(datadir))
    assert os.path.exists("{0}/outdir/test.ucsc.vcf".format(datadir))

def test_drop_genotypes(vcf_object, datadir, clean_tmp):
    """
    Test the method to drop the genotype information from a VCF file
    """
    vcf_object.drop_genotypes(outfile="{0}/outdir/test.sites.vcf.gz".format(datadir), verbose=True)

    assert os.path.exists("{0}/outdir/test.sites.vcf.gz".format(datadir))

def test_drop_info(vcf_object, datadir, clean_tmp):
    """
    Test the method to drop the INFO annotation from a VCF file
    """
    vcf_object.drop_info(outfile="{0}/outdir/test.noinfo.vcf.gz".format(datadir), verbose=True)

    assert os.path.exists("{0}/outdir/test.noinfo.vcf.gz".format(datadir))

def test_convert_PL2GL(datadir, bcftools_folder, clean_tmp):
    '''
    Test the method change PL fields to GL in a VCF file
    '''
    vcf_object = VcfUtils(vcf="{0}/test.gatk.vcf.gz".format(datadir),
                          bcftools_folder=bcftools_folder)

    vcf_object.convert_PL2GL(outfile="{0}/outdir/test.gatk.GL.vcf.gz".format(datadir), verbose=True)

    assert os.path.exists("{0}/outdir/test.gatk.GL.vcf.gz".format(datadir))
