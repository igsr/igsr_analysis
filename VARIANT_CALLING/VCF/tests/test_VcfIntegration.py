import os
import pytest
import glob
import warnings
from VcfIntegration import VcfIntegration

# test_vcfqc.py
@pytest.fixture
def vcf_object(scope='module'):
    '''Returns an  object'''
    print("Creating the object\n")
    vcf_file = pytest.config.getoption("--vcf")
    snptools_folder = pytest.config.getoption("--snptools_folder")
    beagle_folder = pytest.config.getoption("--beagle_folder")
    vcf_object=VcfIntegration(vcf=vcf_file,snptools_folder=snptools_folder,
                           beagle_folder=beagle_folder)
    return vcf_object

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/SNPTools/outdir/*')
    for f in files:
        os.remove(f)
    files1 = glob.glob('data/BEAGLE/outdir/*')
    for f1 in files1:
        os.remove(f1)


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
                          outdir="data/BEAGLE/outdir/",
                          window=12000,
                          overlap=2000,
                          niterations=15,
                          verbose=True)
    assert os.path.exists("data/BEAGLE/outdir/NA12878_chr1_1000000_1001000.beagle.vcf.gz")

def test_run_beagle_multithreaded(vcf_object):
    vcf_object.run_beagle(outprefix="NA12878_chr1_1000000_1001000_mts",
                          outdir="data/BEAGLE/outdir/",
                          window=12000,
                          overlap=2000,
                          niterations=15,
                          nthreads=2)
    assert os.path.exists("data/BEAGLE/outdir/NA12878_chr1_1000000_1001000_mts.beagle.vcf.gz")

def test_run_beagle_with_region(vcf_object):
    vcf_object.run_beagle(outprefix="NA12878_chr1_1000000_1001000",
                          outdir="data/BEAGLE/outdir/",
                          window=12000,
                          overlap=2000,
                          niterations=15,
                          region="chr1:10000-50000")
    assert os.path.exists("data/BEAGLE/outdir/NA12878_chr1_1000000_1001000.chr1:10000-50000.beagle.vcf.gz")

def test_make_beagle_chunks():
    vcf_object=VcfIntegration(vcf="data/BEAGLE/GLs.HG00136.vcf.gz",
                           makeBGLCHUNKS_folder=pytest.config.getoption("--makeBGLCHUNKS_folder"))
    outfile=vcf_object.make_beagle_chunks(window=100,overlap=20,outfile='data/BEAGLE/outdir/chunks.coords')
    assert os.path.exists(outfile)

def test_make_beagle_chunks_with_ucsc_correction():
    vcf_object=VcfIntegration(vcf="data/BEAGLE/GLs.HG00136.ucsc.vcf.gz",
                              makeBGLCHUNKS_folder=pytest.config.getoption("--makeBGLCHUNKS_folder"))
    outfile=vcf_object.make_beagle_chunks(window=100,overlap=20,outfile='data/BEAGLE/outdir/chunks.coords',
                                          correct=True,chrname='chr22')
    assert os.path.exists(outfile)

def test_prepareGenFromBeagle4():
    vcf_object=VcfIntegration(vcf="data/BEAGLE/GLs.HG00136.vcf.gz",
                              makeBGLCHUNKS_folder=pytest.config.getoption("--makeBGLCHUNKS_folder"))
    outfile=vcf_object.make_beagle_chunks(window=100,overlap=20,outfile='data/BEAGLE/outdir/chunks.coords')
    assert os.path.exists(outfile)
