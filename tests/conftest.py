import pytest
import os
import shutil
import glob

from VariantCalling import GATK

@pytest.fixture(autouse=True)
def datadir():
    """
    Returns directory containing data for tests
    """
    return os.path.abspath('./data/')

@pytest.fixture
def vcflib_folder(datadir):
    """Returns the folder containing the vcflib-related binaries"""

    vcflib_folder = None
    if shutil.which('vcf2fasta') is None:
        raise Exception("'vcf2fasta' needs to by in $PATH")
    else:
        vcflib_folder = os.path.dirname(shutil.which('vcf2fasta'))
    
    return vcflib_folder

@pytest.fixture
def vt_folder(datadir):
    """Returns the folder containing the vt binary"""

    vt_folder = None
    if shutil.which('vt') is None:
        raise Exception("'vt' needs to by in $PATH")
    else:
        vt_folder = os.path.dirname(shutil.which('vt'))
    
    return vt_folder


@pytest.fixture
def gatk_folder(datadir):
    '''Returns the folder containing the gatk wrapper script'''

    gatk_folder = None # folder containing the gatk wrapper script
    if shutil.which('gatk3') is None:
        raise Exception("'gatk3' needs to by in $PATH")
    else:
        gatk_folder = os.path.dirname(shutil.which('gatk3'))
    
    return gatk_folder

@pytest.fixture
def gatk_jar_folder():
    '''Returns folder containin the GATK jar file'''

    return '/nfs/production/reseq-info/work/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/'

@pytest.fixture
def picard_folder():
    '''Returns the folder containing the picard .jar file'''
    
    return '~/bin/'


@pytest.fixture
def bgzip_folder(datadir):
    '''Returns the folder containing the bgzip binary'''

    bgzip_folder = None
    if shutil.which('bgzip') is None:
        raise Exception("'bgzip' needs to by in $PATH")
    else:
        bgzip_folder = os.path.dirname(shutil.which('bgzip'))
    
    return bgzip_folder

@pytest.fixture
def gatk_object(gatk_folder, datadir, bgzip_folder):
    '''Returns a GATK object'''

    bam_file = "{0}/exampleBAM.bam".format(datadir)
    reference = "{0}/exampleFASTA.fasta".format(datadir)
    
    gatk_object = GATK(bam=bam_file, reference=reference,
                       gatk_folder=gatk_folder, bgzip_folder=bgzip_folder)

    return gatk_object

@pytest.fixture
def clean_tmp(datadir):
    yield
    print("Cleanup files")
    files = glob.glob("{0}/outdir/*".format(datadir))
    for f in files:
        os.remove(f)

@pytest.fixture
def bcftools_folder(datadir):
    '''Returns the folder containing the bcftools binary'''

    bcftools_folder = None
    if shutil.which('bcftools') is None:
        raise Exception("'bcftools needs to by in $PATH")
    else:
        bcftools_folder = os.path.dirname(shutil.which('bcftools'))
    
    return bcftools_folder

@pytest.fixture
def bedtools_folder(datadir):
    '''Returns the folder containing the bedtools binary'''

    bedtools_folder = None
    if shutil.which('bedtools') is None:
        raise Exception("'bedtools needs to by in $PATH")
    else:
        bedtools_folder = os.path.dirname(shutil.which('bedtools'))
    
    return bedtools_folder

@pytest.fixture
def shapeit_folder(datadir):
    '''Returns the folder containing the shapeit binary'''

    shapeit_folder = None
    if shutil.which('shapeit') is None:
        raise Exception("'shapeit' needs to by in $PATH")
    else:
        shapeit_folder = os.path.dirname(shutil.which('shapeit'))
    
    return shapeit_folder

@pytest.fixture
def hive_dir():
    '''Returns directory containing Ensembl ehive codebase'''
    
    assert os.getenv('ENSEMBL_HIVE_DIR'), "$ENSEMBL_HIVE_DIR undefined. You need to set this ENV variable to run this test suite"
    return os.getenv('ENSEMBL_HIVE_DIR')

@pytest.fixture
def beagle_bins():
    '''Returns a tuple with directory containing the Beagle .jar file and Beagle's .jar filename'''

    return '~/bin/beagle/', 'beagle.08Jun17.d8b.jar'
