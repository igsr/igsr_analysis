import os
import glob
import pytest

from VCF.VCFIntegration import Beagle

# test_VcfIntegration_Beagle.py
@pytest.fixture
def vcf_object(scope='module'):
    '''Returns an  object'''
    print("Creating the object\n")
    vcf_file = pytest.config.getoption("--vcf")
    beagle_folder = pytest.config.getoption("--beagle_folder")
    beagle_jar = pytest.config.getoption("--beagle_jar")
    vcf_object = Beagle(vcf=vcf_file, beagle_jar=beagle_jar,
                        beagle_folder=beagle_folder)
    return vcf_object

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/outdir/*')
    for f in files:
        os.remove(f)

def test_run_beagle(vcf_object):
    outfile = vcf_object.run_beagle(outprefix="NA12878_chr1_1000000_1001000",
                                    outdir="data/outdir/",
                                    window=12000,
                                    overlap=2000,
                                    niterations=15,
                                    verbose=True)
    assert os.path.exists(outfile)

def test_run_beagle_multithreaded(vcf_object):
    vcf_object.run_beagle(outprefix="NA12878_chr1_1000000_1001000_mts",
                          outdir="data/outdir/",
                          window=12000,
                          overlap=2000,
                          niterations=15,
                          nthreads=2)
    assert os.path.exists("data/outdir/NA12878_chr1_1000000_1001000_mts.beagle.vcf.gz")

def test_make_beagle_chunks():
    vcf_object = Beagle(vcf="data/BEAGLE/GLs.HG00136.vcf.gz",
                        makeBGLCHUNKS_folder=pytest.config.getoption("--makeBGLCHUNKS_folder"))

    outfile = vcf_object.make_beagle_chunks(window=700,
                                            overlap=200,
                                            outfile='data/outdir/chunks.coords',
                                            verbose=True)
    assert os.path.exists(outfile)

def test_run_beagle_with_region_correct():
    '''
    Run Beagle on all chunks created in the previous test
    '''
    vcf_object = Beagle(vcf="data/BEAGLE/GLs.HG00136.vcf.gz",
                        beagle_folder=pytest.config.getoption("--beagle_folder"),
                        beagle_jar=pytest.config.getoption("--beagle_jar"))

    with open('data/outdir/chunks.coords') as f:
        for line in f:
            line = line.rstrip('\n')
            chroname = line.split('\t')[0]
            start = line.split('\t')[1]
            end = line.split('\t')[2]
            beagle_out = vcf_object.run_beagle(outprefix="GLs.HG00136.correct",
                                               region="{0}:{1}-{2}".format(chroname, start, end),
                                               outdir="data/outdir/",
                                               correct=True,
                                               verbose=True)
    assert os.path.exists(beagle_out)

def test_prepareGenFromBeagle4(clean_tmp):
    vcf_object = Beagle(vcf="data/GLs.HG00136.vcf.gz",
                        prepareGenFromBeagle4_folder=pytest.config.getoption("--prepareGenFromBeagle4_folder"))

    outdict = vcf_object.prepare_Gen_From_Beagle4(prefix_in="data/outdir/GLs.HG00136.correct.22",
                                                  outprefix='data/outdir/input.shapeit.22',
                                                  verbose=True)
    assert os.path.exists(outdict['gen_gz'])
    assert os.path.exists(outdict['gen_sample'])
    assert os.path.exists(outdict['hap_gz'])
    assert os.path.exists(outdict['hap_sample'])
