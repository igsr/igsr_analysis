import os
import pytest
import glob
from VCF.VCFIntegration import SNPTools

# test_vcfqc.py
@pytest.fixture
def vcf_object(scope='module'):
    '''Returns an  object'''
    print("Creating the object\n")
    vcf_file = pytest.config.getoption("--vcf")
    snptools_folder = pytest.config.getoption("--snptools_folder")
    vcf_object = SNPTools(vcf=vcf_file, snptools_folder=snptools_folder)
    return vcf_object

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/SNPTools/outdir/*')
    for f in files:
        os.remove(f)

def test_run_bamodel(vcf_object):
    vcf_object.run_bamodel(sample="NA12878_chr1_1000000_1001000",
                           bamfiles="data/SNPTools/bamlist.txt",
                           outdir="data/SNPTools/outdir/")
    assert os.path.exists("data/SNPTools/outdir/NA12878_chr1_1000000_1001000.raw")

def test_run_poprob(vcf_object):
    vcf_object.run_poprob(outprefix="NA12878_chr1_1000000_1001000",
                          rawlist="data/SNPTools/rawlist.txt",
                          outdir="data/SNPTools/outdir/")
    assert os.path.exists("data/SNPTools/outdir/NA12878_chr1_1000000_1001000.prob")

def test_run_prob2vcf(vcf_object):
    vcf_object.run_prob2vcf(probf="data/SNPTools/outdir/NA12878_chr1_1000000_1001000.prob",
                            outprefix="NA12878_chr1_1000000_1001000",
                            outdir="data/SNPTools/outdir/",
                            chro="chr1")
    assert os.path.exists("data/SNPTools/outdir/NA12878_chr1_1000000_1001000.vcf.gz")
