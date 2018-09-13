import os
import pytest
import subprocess
import glob
from PyHive.VcfIntegration import *

# test_pyhive_runBeagle.py

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/outdir/*')
    for f in files:
        os.remove(f)

def test_runBeagle():
    vcf_f= pytest.config.getoption("vcf_chr20")
    beagle_folder=pytest.config.getoption("beagle_folder")
    beagle_jar=pytest.config.getoption("beagle_jar")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/" 
    work_dir= "data/outdir"

    command="perl {0}/standaloneJob.pl PyHive.VcfIntegration.run_Beagle -language python3 -vcf_file {1} -beagle_folder {2} -beagle_jar {3} -work_dir {4} -outprefix {5} \
    -window {6} -overlap {7} -niterations {8} -nthreads {9} -region_chunk {10} -verbose True".format(hive_scripts, vcf_f, beagle_folder, beagle_jar, work_dir, 
                                                                                                    'test', 12000, 2000, 15, 1, 'chr20:1000000-1000150')
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_runBeagle_noparams():
    vcf_f= pytest.config.getoption("vcf_chr20")
    beagle_folder=pytest.config.getoption("beagle_folder")
    beagle_jar=pytest.config.getoption("beagle_jar")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"
    work_dir= "data/outdir"

    command="perl {0}/standaloneJob.pl PyHive.VcfIntegration.run_Beagle -language python3 -vcf_file {1} -beagle_folder {2} -beagle_jar {3} -work_dir {4} -outprefix {5} \
    -region_chunk {6} -verbose True".format(hive_scripts, vcf_f, beagle_folder, beagle_jar, work_dir,
                                        'test', 'chr20:1000000-1000150')

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_runBeagle_noparams_correct(clean_tmp):
    vcf_f= pytest.config.getoption("vcf_chr20")
    beagle_folder=pytest.config.getoption("beagle_folder")
    beagle_jar=pytest.config.getoption("beagle_jar")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"
    work_dir= "data/outdir"

    command="perl {0}/standaloneJob.pl PyHive.VcfIntegration.run_Beagle -language python3 -vcf_file {1} -beagle_folder {2} -beagle_jar {3} -work_dir {4} -outprefix {5} \
    -region_chunk {6} -correct True -verbose True".format(hive_scripts, vcf_f, beagle_folder, beagle_jar, work_dir,
                                                          'test_correct', 'chr20:1000000-1000150')
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)



