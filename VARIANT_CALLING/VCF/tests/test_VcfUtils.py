import os
import pytest
import glob
from VcfUtils import VcfUtils

# test_VcfUtils.py

@pytest.fixture
def vcf_object():
    '''Returns an  object'''
    vcf_file = pytest.config.getoption("--vcf")
    vcflist = pytest.config.getoption("--vcflist")
    bcftools_folder = pytest.config.getoption("--bcftools_folder")
    gatk_folder = pytest.config.getoption("--gatk_folder")
    vcf_object=VcfUtils(vcf=vcf_file,vcflist=vcflist,bcftools_folder=bcftools_folder,
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
    files = glob.glob('data/out/*')
    for f in files:
        os.remove(f)

def test_correct_ambiguity(vcf_ambiguity):
    outfile=vcf_ambiguity.correct_ambiguity_codes(outfile='data/out/test.corrected.vcf.gz')

    assert os.path.exists("data/out/test.corrected.vcf.gz")

def test_vcf_reheader(vcf_object):
    outfile=vcf_object.reheader(newheader="data/newheader.txt", outprefix="data/out/test1") 

    assert os.path.exists("data/out/test1.reheaded.vcf.gz")

def test_vcf_reheader_with_samplef(vcf_object):
    '''
    Test the reheader method and add new sample names
    '''
    outfile=vcf_object.reheader(newheader="data/newheader.txt", samplefile="data/samples.txt",
                                outprefix="data/out/test2")

    assert os.path.exists("data/out/test2.reheaded.vcf.gz")

def test_combine_uncompressed(vcf_object):
    '''
    Test the combine method producing a VCF
    '''
    vcf_object.combine(labels=['gatk','lc_bcftools'],reference='data/exampleFASTA.fasta',outprefix='out_combine',
                       outdir='data/out/',genotypemergeoption='UNIQUIFY')

    assert os.path.exists("data/out/out_combine.vcf")
    assert os.path.exists("data/out/out_combine.vcf.idx")

def test_combine_compressed(vcf_object):
    '''
    Test the combine method producing a VCF.gz file and passing also 
    some options
    '''
    vcf_object.combine(labels=['gatk','lc_bcftools'],reference='data/exampleFASTA.fasta',outprefix='out_combine',
                       outdir='data/out/',compress=True,genotypemergeoption='UNIQUIFY',options=['-env','-sites_only',
                                                                                                '--filteredAreUncalled'])
    
    assert os.path.exists("data/out/out_combine.vcf.gz")

def test_change_chrnames_2ensembl(vcf_object):
    '''
    Test the method to change the style of the chrnames (from UCSC to Ensembl)
    '''
    vcf_object.rename_chros(chr_types='ensembl', outfile='data/out/test.ensembl.vcf.gz')

    assert os.path.exists("data/out/test.ensembl.vcf.gz")

def test_change_chrnames_2ucsc(clean_tmp):
    '''
    Test the method to change the style of the chrnames (from Ensembl to UCSC)
    '''
    vcf_object=VcfUtils(vcf='data/out/test.ensembl.vcf.gz')

    vcf_object.rename_chros(chr_types='ucsc', outfile='data/out/test.ucsc.vcf.gz')

    assert os.path.exists("data/out/test.ucsc.vcf.gz")

