import os
import pytest
import subprocess
import glob

from VcfUtils import VcfUtils

# test_pyhive_Vcf.py

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('out/*')
    for f in files:
        os.remove(f)

def test_VcfReheader():
    '''
    Test PyHive.Vcf.VcfReheader after setting the 'work_dir' arg
    '''
    vcf_file = pytest.config.getoption("--vcf")
    bcftools_folder = pytest.config.getoption("--bcftools_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"

    command="perl {0}/standaloneJob.pl PyHive.Vcf.VcfReheader -language python3 -filepath {1} -bcftools_folder {2} -newheader data/newheader.txt -work_dir data/outdir/ -debug 1".format(hive_scripts, vcf_file, bcftools_folder)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_VcfReheader_w_samplename():
    '''
    Test PyHive.Vcf.VcfReheader after setting the 'samplename' arg
    '''
    vcf_file = pytest.config.getoption("--vcf")
    bcftools_folder = pytest.config.getoption("--bcftools_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"

    command="perl {0}/standaloneJob.pl PyHive.Vcf.VcfReheader -language python3 -filepath {1} -bcftools_folder {2} -newheader data/newheader.txt -work_dir data/outdir/ -samplename testname -debug 1".format(hive_scripts, vcf_file, bcftools_folder)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_VcfReheader_w_samplenfile():
    '''
    Test PyHive.Vcf.VcfReheader after setting the 'samplename' arg
    '''
    vcf_file = pytest.config.getoption("--vcf")
    bcftools_folder = pytest.config.getoption("--bcftools_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"

    command="perl {0}/standaloneJob.pl PyHive.Vcf.VcfReheader -language python3 -filepath {1} -bcftools_folder {2} -newheader data/newheader.txt -work_dir data/outdir/ -samplefile data/samplefile.txt -debug 1".format(hive_scripts, vcf_file, bcftools_folder)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_VcfReplaceChrNames():
    '''
    Test PyHive.Vcf.VcfReplaceChrNames
    '''
    vcf_file = pytest.config.getoption("--vcf")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"
    bgzip_folder= pytest.config.getoption("--bgzip_folder")

    command="perl {0}/standaloneJob.pl PyHive.Vcf.VcfReplaceChrNames -language python3 -filepath {1} -chr_types 'ensembl' -work_dir data/outdir/ -outprefix test.ensembl -bgzip_folder {2}".format(hive_scripts, vcf_file, bgzip_folder)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_VcfconvertPL2GL(clean_tmp):
    '''
    Test PyHive.Vcf.convertPL2GL
    '''
    vcf_file = 'data/test.gatk.vcf.gz'
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"
    bcftools_folder= pytest.config.getoption("--bcftools_folder")

    command="perl {0}/standaloneJob.pl PyHive.Vcf.convertPL2GL -language python3 -filepath {1} -work_dir data/outdir/ -outprefix test.pl2gl -bcftools_folder {2}".format(hive_scripts, vcf_file, bcftools_folder)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)
