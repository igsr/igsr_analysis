import os
import pytest
import subprocess
from VcfUtils import VcfUtils

# test_Vcf.py

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    os.remove('data/test.vcf.gz.reheaded.vcf.gz')


def test_VcfReheader(clean_tmp):
    '''
    Test PyHive.Vcf.VcfReheader with the minimum number of args that 
    are possible
    '''
    vcf_file = pytest.config.getoption("--vcf")
    bcftools_folder = pytest.config.getoption("--bcftools_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/" 

    command="perl {0}/standaloneJob.pl PyHive.Vcf.VcfReheader -language python3 -filepath {1} -bcftools_folder {2} -newheader data/newheader.txt -debug 1".format(hive_scripts, vcf_file, bcftools_folder)
    
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)
    
def test_VcfReheader_w_workingdir(clean_tmp):
    '''
    Test PyHive.Vcf.VcfReheader after setting the 'work_dir' arg
    '''
    vcf_file = pytest.config.getoption("--vcf")
    bcftools_folder = pytest.config.getoption("--bcftools_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"

    command="perl {0}/standaloneJob.pl PyHive.Vcf.VcfReheader -language python3 -filepath {1} -bcftools_folder {2} -newheader data/newheader.txt -work_dir data/ -debug 1".format(hive_scripts, vcf_file, bcftools_folder)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_VcfReheader_w_samplename(clean_tmp):
    '''
    Test PyHive.Vcf.VcfReheader after setting the 'samplename' arg
    '''
    vcf_file = pytest.config.getoption("--vcf")
    bcftools_folder = pytest.config.getoption("--bcftools_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"

    command="perl {0}/standaloneJob.pl PyHive.Vcf.VcfReheader -language python3 -filepath {1} -bcftools_folder {2} -newheader data/newheader.txt -work_dir data -samplename testname -debug 1".format(hive_scripts, vcf_file, bcftools_folder)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_VcfReheader_w_samplenfile(clean_tmp):
    '''
    Test PyHive.Vcf.VcfReheader after setting the 'samplename' arg
    '''
    vcf_file = pytest.config.getoption("--vcf")
    bcftools_folder = pytest.config.getoption("--bcftools_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"

    command="perl {0}/standaloneJob.pl PyHive.Vcf.VcfReheader -language python3 -filepath {1} -bcftools_folder {2} -newheader data/newheader.txt -work_dir data -samplefile data/samplefile.txt -debug 1".format(hive_scripts, vcf_file, bcftools_folder)
    print(command)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)
