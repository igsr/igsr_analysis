import os
import pytest
import subprocess
import glob
from VcfGenotype import VcfGenotype

# test_runBeagle.py

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/outdir/*')
    for f in files:
        os.remove(f)

def test_runBeagle(clean_tmp):
    vcf_f= pytest.config.getoption("vcf")
    beagle_folder=pytest.config.getoption("beagle_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/" 
    work_dir= "data/outdir"

    command="perl {0}/standaloneJob.pl PyHive.VcfGenotype.run_Beagle -language python3 -vcf_file {1} -beagle_folder {2} -work_dir {3} -outprefix {4} \
             -window {5} -overlap {6} -niterations {7} -verbose True".format(hive_scripts, vcf_f, beagle_folder, work_dir, 'test', 12000, 2000, 15 )
    
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)
    
