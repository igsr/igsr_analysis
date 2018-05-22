import os
import pytest
import glob
import pdb

from VcfUtils import VcfUtils

# test_VcfUtils.py

@pytest.fixture
def vcf_object():
    '''Returns an  object'''

    vcf_file = pytest.config.getoption("--vcf")
    vcflist = pytest.config.getoption("--vcflist")
    bcftools_folder = pytest.config.getoption("--bcftools_folder")
    bgzip_folder = pytest.config.getoption("--bgzip_folder")
    gatk_folder = pytest.config.getoption("--gatk_folder")
    vcf_object=VcfUtils(vcf=vcf_file,bgzip_folder=bgzip_folder,vcflist=vcflist,bcftools_folder=bcftools_folder,
                        gatk_folder=gatk_folder)
    return vcf_object

@pytest.fixture
def vcf_ambiguity():
    '''
    Returns a VcfUtils object that contains the REF or ALT column with some ambiguity codes
    '''
    vcf_file = pytest.config.getoption("--vcf_ambiguity")
    vcf_object=VcfUtils(vcf=vcf_file)
    return vcf_object


@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/outdir/*')
    for f in files:
        os.remove(f)

def test_correct_ambiguity(vcf_ambiguity):

    outfile=vcf_ambiguity.correct_ambiguity_codes(outfile='data/outdir/test.corrected.vcf.gz')

    assert os.path.exists("data/outdir/test.corrected.vcf.gz")
def test_vcf_reheader(vcf_object):
    outfile=vcf_object.reheader(newheader="data/newheader.txt", outprefix="data/outdir/test1", verbose=True) 

    assert os.path.exists("data/outdir/test1.reheaded.vcf.gz")


def test_vcf_reheader_with_samplef(vcf_object):
    '''
    Test the reheader method and add new sample names
    '''
    outfile=vcf_object.reheader(newheader="data/newheader.txt", samplefile="data/samples.txt",
                                outprefix="data/outdir/test2")

    assert os.path.exists("data/outdir/test2.reheaded.vcf.gz")

def test_combine_uncompressed(vcf_object):
    '''
    Test the combine method producing a VCF
    '''
    vcf_object.combine(labels=['gatk','lc_bcftools'],reference='data/exampleFASTA.fasta',outprefix='out_combine',
                       outdir='data/outdir/', verbose=True, genotypemergeoption='UNIQUIFY')

    assert os.path.exists("data/outdir/out_combine.vcf")
    assert os.path.exists("data/outdir/out_combine.vcf.idx")

def test_combine_compressed(vcf_object):
    '''
    Test the combine method producing a VCF.gz file and passing also 
    some options
    '''
    vcf_object.combine(labels=['gatk','lc_bcftools'],reference='data/exampleFASTA.fasta',outprefix='out_combine',
                       outdir='data/outdir/',compress=True,genotypemergeoption='UNIQUIFY',options=['-env','-sites_only',
                                                                                           '--filteredAreUncalled'])
    
    assert os.path.exists("data/outdir/out_combine.vcf.gz")

def test_change_chrnames_2ensembl(vcf_object):
    '''
    Test the method to change the style of the chrnames (from UCSC to Ensembl)
    '''

    vcf_object.rename_chros(chr_types='ensembl', outfile='data/outdir/test.ensembl.vcf.gz')
    vcf_object.rename_chros(chr_types='ensembl', compress=False, outfile='data/outdir/test.ensembl.vcf')

    assert os.path.exists("data/outdir/test.ensembl.vcf.gz")
    assert os.path.exists("data/outdir/test.ensembl.vcf")

def test_change_chrnames_2ucsc():
    '''
    Test the method to change the style of the chrnames (from Ensembl to UCSC)
    '''
    vcf_object=VcfUtils(vcf='data/outdir/test.ensembl.vcf.gz',
                        bgzip_folder=pytest.config.getoption("--bgzip_folder"))

    vcf_object.rename_chros(chr_types='ucsc', outfile='data/outdir/test.ucsc.vcf.gz')
    vcf_object.rename_chros(chr_types='ucsc', outfile='data/outdir/test.ucsc.vcf', compress=False)
    
    assert os.path.exists("data/outdir/test.ucsc.vcf.gz")
    assert os.path.exists("data/outdir/test.ucsc.vcf")

def test_drop_genotypes(vcf_object):
    '''
    Test the method to drop the genotype information from a VCF file
    '''
    vcf_object.drop_genotypes(outfile='data/outdir/test.sites.vcf.gz',verbose=True)

    assert os.path.exists("data/outdir/test.sites.vcf.gz")

def test_drop_info(vcf_object):
    '''
    Test the method to drop the INFO annotation from a VCF file
    '''
    vcf_object.drop_info(outfile='data/outdir/test.noinfo.vcf.gz',verbose=True)

    assert os.path.exists("data/outdir/test.noinfo.vcf.gz")

def test_convert_PL2GL(clean_tmp):
    '''
    Test the method change PL fields to GL in a VCF file
    '''
    vcf_object=VcfUtils(vcf='data/test.gatk.vcf.gz',
                        bcftools_folder=pytest.config.getoption("--bcftools_folder"))
    
    vcf_object.convert_PL2GL(outfile='data/outdir/test.gatk.GL.vcf.gz',verbose=True)

    assert os.path.exists("data/outdir/test.gatk.GL.vcf.gz")

