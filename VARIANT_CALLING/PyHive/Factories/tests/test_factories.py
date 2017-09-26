import os
import pytest
import subprocess
import glob

# test_factories.py

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/outdir/*')
    for f in files:
        os.remove(f)

def test_chrfactory():
    fa_ix= pytest.config.getoption("faix")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/" 

    command="perl {0}/standaloneJob.pl PyHive.Factories.ChrFactory -language python3 -faix {1}".format(hive_scripts, fa_ix)
    
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)    

def test_BeagleChunkFactory(clean_tmp):
    vcf_f= pytest.config.getoption("vcf_gts")
    beaglechunks_folder=pytest.config.getoption("makeBGLCHUNKS_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"
    work_dir= "data/outdir"

    command="perl {0}/standaloneJob.pl PyHive.Factories.BeagleChunkFactory -language python3 -filepath {1} -makeBGLCHUNKS_folder {2} -work_dir {3} -window {4} \
    -overlap {5}".format(hive_scripts, vcf_f, beaglechunks_folder, work_dir, 100, 2)
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)


