import os
import glob
import pytest
import shutil
import pdb

from VCF.VCFIntegration import Beagle

# test_VcfIntegration_Beagle.py
@pytest.fixture
def beagle_object(beagle_bins, datadir):
    """Returns a Beagle object"""

    vcf_file = "{0}/test.vcf.gz".format(datadir)
    beagle_folder = beagle_bins[0]
    beagle_jar = beagle_bins[1]

    beagleObj = Beagle(vcf=vcf_file, beagle_jar=beagle_jar,
                       beagle_folder=beagle_folder)
    
    return beagleObj

def test_run_beagle(beagle_object, datadir, clean_tmp):

    outfile = beagle_object.run_beagle(outprefix="NA12878_chr1_1000000_1001000",
                                       outdir="{0}/outdir/".format(datadir),
                                       window=12000,
                                       overlap=2000,
                                       niterations=15,
                                       verbose=True)

    assert os.path.exists(outfile)

def test_run_beagle_multithreaded(beagle_object, datadir, clean_tmp):

    beagle_object.run_beagle(outprefix="NA12878_chr1_1000000_1001000_mts",
                          outdir="{0}/outdir/".format(datadir),
                          window=12000,
                          overlap=2000,
                          niterations=15,
                          nthreads=2)

    assert os.path.exists("{0}/outdir/NA12878_chr1_1000000_1001000_mts.beagle.vcf.gz".format(datadir))

def test_make_beagle_chunks(datadir):

    makeBGLCHUNKS_folder = None # folder containing the makeBGLCHUNKS binary
    if shutil.which('makeBGLCHUNKS') is None:
        raise Exception("'makeBGLCHUNKS' needs to by in $PATH")
    else:
        makeBGLCHUNKS_folder = os.path.dirname(shutil.which('makeBGLCHUNKS'))
    
    vcf_object = Beagle(vcf="{0}/BEAGLE/GLs.HG00136.vcf.gz".format(datadir),
                        makeBGLCHUNKS_folder=makeBGLCHUNKS_folder)

    outfile = vcf_object.make_beagle_chunks(window=700,
                                            overlap=200,
                                            outfile="{0}/outdir/chunks.coords".format(datadir),
                                            verbose=True)
    assert os.path.exists(outfile)

def test_run_beagle_with_region_correct(datadir, beagle_bins, clean_tmp):
    '''
    Run Beagle on all chunks created in the previous test
    '''
    vcf_object = Beagle(vcf="{0}/BEAGLE/GLs.HG00136.vcf.gz".format(datadir),
                        beagle_folder=beagle_bins[0],
                        beagle_jar=beagle_bins[1])

    with open("{0}/outdir/chunks.coords".format(datadir)) as f:
        for line in f:
            line = line.rstrip('\n')
            chroname = line.split('\t')[0]
            start = line.split('\t')[1]
            end = line.split('\t')[2]
            beagle_out = vcf_object.run_beagle(outprefix="GLs.HG00136.correct",
                                               region="{0}:{1}-{2}".format(chroname, start, end),
                                               outdir="{0}/outdir/".format(datadir),
                                               correct=True,
                                               verbose=True)
    assert os.path.exists(beagle_out)

def test_prepareGenFromBeagle4(datadir, clean_tmp):
 
    prepareGenFromBeagle4_folder = None # folder containing the prepareGenFromBeagle4 binary
    if shutil.which('prepareGenFromBeagle4') is None:
        raise Exception("'prepareGenFromBeagle4' needs to by in $PATH")
    else:
        prepareGenFromBeagle4_folder = os.path.dirname(shutil.which('prepareGenFromBeagle4'))

    vcf_object = Beagle(vcf="{0}/GLs.HG00136.vcf.gz".format(datadir),
                        prepareGenFromBeagle4_folder=prepareGenFromBeagle4_folder)

    outdict = vcf_object.prepare_Gen_From_Beagle4(prefix_in="{0}/outdir/GLs.HG00136.correct.22".format(datadir),
                                                  outprefix="{0}/outdir/input.shapeit.22".format(datadir),
                                                  verbose=True)
    assert os.path.exists(outdict['gen_gz'])
    assert os.path.exists(outdict['gen_sample'])
    assert os.path.exists(outdict['hap_gz'])
    assert os.path.exists(outdict['hap_sample'])
