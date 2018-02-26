import os
import pytest
import subprocess
import glob

from GATK import GATK

# test_pyhive_runGATK_UG.py

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('out/*')
    for f in files:
        os.remove(f)

def test_runGATK_UG(clean_tmp):

    shapeit_folder=pytest.config.getoption("shapeit_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"
    bam_file = pytest.config.getoption("--bam")
    reference = pytest.config.getoption("--reference")
    gatk_folder= pytest.config.getoption("--gatk_folder")
    bgzip_folder = pytest.config.getoption("--bgzip_folder")
    glm = pytest.config.getoption("--glm")
    output_mode = pytest.config.getoption("--output_mode")

    work_dir= "out/"

    command="perl {0}/standaloneJob.pl PyHive.VariantCalling.GATK_UG -language python3 \
    -outprefix {1} -work_dir {2} -chunk {3} -bamlist {4} -reference {5} \
    -gatk_folder {6} -bgzip_folder {7} -glm {8} -output_mode {9} -verbose True".format(hive_scripts, 'out', work_dir, 
                                                                                       'chr1:10000-30000', bam_file, 
                                                                                       reference, gatk_folder, bgzip_folder,
                                                                                       glm, output_mode)
    
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

