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
def gatk_folder(datadir):
    '''Returns the folder containing the gatk wrapper script'''

    gatk_folder = None # folder containing the gatk wrapper script
    if shutil.which('gatk3') is None:
        raise Exception("'gatk3' needs to by in $PATH")
    else:
        gatk_folder = os.path.dirname(shutil.which('gatk3'))
    
    return gatk_folder

@pytest.fixture
def bgzip_folder(datadir):
    '''Returns the folder containing the bgzip binary'''

    bgzip_folder = None # folder containing bgzip
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

    shapeit_folder = None # folder containing the Shapeit binary
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


"""
def pytest_addoption(parser):
    parser.addoption('--hive_lib', default='~/lib/ensembl-hive_2.4/' ,action='store_true', help='Path folder containing eHive scripts')
    parser.addoption('--vcf', default='data/test.vcf.gz' ,action='store_true', help='Path to vcf file')
    parser.addoption('--region', default='data/region.bed' ,action='store_true', help='BED file with a small region in chr1')
    parser.addoption('--vcf_gts_ucsc', default='data/GLs.HG00136.ucsc.vcf.gz' ,action='store_true', help='Path to vcf file with GTs with UCSC-style chro names')
    parser.addoption('--vt_folder', default='~/bin/vt/' ,action='store_true', help='Path to folder containing vt binary')
    parser.addoption('--chk_indel_folder', default='~/bin/' ,action='store_true', help='Folder with chk_indel_rg binary')
    parser.addoption('--samtools_folder', default='/homes/ernesto/bin/samtools-1.6/bin/', action='store_true', help='Folder with samtools binary')
    parser.addoption('--java_folder', default='/usr/bin/',action='store_true', help='Folder with java binary')
    parser.addoption('--picard_folder', default='~/bin/',action='store_true', help='Folder with Picard jar file')
    parser.addoption('--bedtools_folder', default='/homes/ernesto/bin/bedtools-2.25.0/bin/', action='store_true', help='Folder with bedtools binary')
    parser.addoption('--makeBGLCHUNKS_folder', default='~/bin/shapeit2_v2_12/bin/makeBGLCHUNKS/bin/' ,action='store_true', help='Folder with makeBGLCHUNKS binary')
    parser.addoption('--bcftools_folder', default='~/bin/bcftools-1.6/' ,action='store_true', help='Folder containing bcftools binaries')
    parser.addoption('--beagle_jar', default='beagle.08Jun17.d8b.jar' ,action='store_true', help='Name of Beagle jar file')
    parser.addoption('--beagle_folder', default='~/bin/beagle/' ,action='store_true', help='Folder with Beagle jar file')
    parser.addoption('--prepareGenFromBeagle4_folder', default='/homes/ernesto/bin/shapeit2_v2_12/bin/prepareGenFromBeagle4/bin/', action='store_true', help='Folder with prepareGenFromBeagle4 binary')
    parser.addoption('--ligateHAPLOTYPES_folder', default='/homes/ernesto/bin/shapeit2_v2_12/bin/ligateHAPLOTYPES/bin/', action='store_true', help='Folder with ligateHAPLOTYPES binary')
    parser.addoption('--shapeit_folder', default='~/bin/shapeit2_v2_12/bin/' ,action='store_true', help='Folder with SHAPEIT binary')
    parser.addoption('--gatk_folder', default='~/bin/GATK/', action='store_true', help='Path to folder containing the GATK jar file')
    parser.addoption('--bgzip_folder', default='/nfs/production/reseq-info/work/ernesto/bin/anaconda3/bin/', action='store_true', help='Path to folder containing the Bgzip binary')
    parser.addoption('--vcflib_folder', default='~/bin/vcflib/bin/', action='store_true', help='Path to folder containing the vcflib binaries')
    parser.addoption('--hostname', default='mysql-g1kdcc-public', action='store_true', help='host name')
    parser.addoption('--username', default='g1krw', action='store_true', help='user name')
    parser.addoption('--port', default=4197, action='store_true', help='port')
    parser.addoption('--pwd', default='test', action='store_true', help='password')
    parser.addoption('--db', default='g1k_archive_staging_track', action='store_true', help='database name')
    parser.addoption('--reference', default='data/exampleFASTA.fasta', action='store_true', help='Path to Fasta file')
    parser.addoption('--vcf_ambiguity', default='data/test.amb.vcf.gz', action='store_true', help='Path to VCF file containing the REF or ALT column with some ambiguity codes')
    parser.addoption('--vcflist', default=['data/test.vcf.gz','data/test1.vcf.gz'], action='store_true', help='List with VCF paths')
    parser.addoption('--snptools_folder', default='~/bin/snptools/', action='store_true', help='Folder with SNPTools binaries')
    parser.addoption('--chr_file', default='data/chr_file.txt', action='store_true', help='File with chros for get_chros function')
    parser.addoption('--glm', default='SNP', action='store_true', help='--glm option for GATK UG') 
    parser.addoption('--output_mode', default='EMIT_ALL_SITES', action='store_true', help='--output_mode option for GATK UG')
 """   
