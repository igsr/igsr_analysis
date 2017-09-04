import os
import pytest
import glob
import warnings
from VcfGenotype import VcfGenotype

# test_vcfqc.py
@pytest.fixture
def vcf_object(scope='module'):
    '''Returns an  object'''
    print("Creating the object\n")
    vcf_file = pytest.config.getoption("--vcf")
    snptools_folder = pytest.config.getoption("--snptools_folder")
    beagle_folder = pytest.config.getoption("--beagle_folder")
    vcf_object=VcfGenotype(vcf=vcf_file,snptools_folder=snptools_folder,
                           beagle_folder=beagle_folder)
    return vcf_object

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/SNPTools/outdir/*')
    for f in files:
        os.remove(f)



def test_run_snptools_bamodel(vcf_object):
    vcf_object.run_snptools_bamodel(sample="NA12878_chr1_1000000_1001000", 
                                    bamfiles="data/SNPTools/bamlist.txt",
                                    outdir="data/SNPTools/outdir/")
    assert os.path.exists("data/SNPTools/outdir/NA12878_chr1_1000000_1001000.raw")

def test_run_snptools_poprob(vcf_object):
    vcf_object.run_snptools_poprob(outprefix="NA12878_chr1_1000000_1001000",
                                   rawlist="data/SNPTools/rawlist.txt",
                                   outdir="data/SNPTools/outdir/")
    assert os.path.exists("data/SNPTools/outdir/NA12878_chr1_1000000_1001000.prob")

def test_run_snptools_prob2vcf(vcf_object):
    vcf_object.run_snptools_prob2vcf(probf="data/SNPTools/outdir/NA12878_chr1_1000000_1001000.prob",
                                     outprefix="NA12878_chr1_1000000_1001000",
                                     outdir="data/SNPTools/outdir/",
                                     chro="chr1")
    assert os.path.exists("data/SNPTools/outdir/NA12878_chr1_1000000_1001000.vcf.gz")

def test_run_beagle(vcf_object):
    vcf_object.run_beagle(outprefix="NA12878_chr1_1000000_1001000",
                          outdir="data/SNPTools/outdir/",
                          window=12000,
                          overlap=2000,
                          niterations=15,
                          verbose=True)
    assert os.path.exists("data/SNPTools/outdir/NA12878_chr1_1000000_1001000.beagle.vcf.gz")

def test_run_beagle_multithreaded(vcf_object, clean_tmp):
    vcf_object.run_beagle(outprefix="NA12878_chr1_1000000_1001000_mts",
                          outdir="data/SNPTools/outdir/",
                          window=12000,
                          overlap=2000,
                          niterations=15,
                          nthreads=2)
    assert os.path.exists("data/SNPTools/outdir/NA12878_chr1_1000000_1001000_mts.beagle.vcf.gz")
