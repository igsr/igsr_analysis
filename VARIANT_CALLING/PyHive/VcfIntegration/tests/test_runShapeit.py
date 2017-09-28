import os
import pytest
import subprocess
import glob
from PyHive.VcfIntegration import *

# test_runShapeit.py

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/outdir/*')
    for f in files:
        os.remove(f)

def test_runShapeit():
    shapeit_folder=pytest.config.getoption("shapeit_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/" 
    work_dir= "data/outdir"

    command="perl {0}/standaloneJob.pl PyHive.VcfIntegration.run_Shapeit -language python3 \
    -shapeit_folder {1} -work_dir {2} -input_gen '{3}' -input_init '{4}' -input_scaffold '{5}' -outprefix {6} \
    -verbose True".format(hive_scripts, shapeit_folder, work_dir, 'data/SHAPEIT/input.shapeit.22.gen.gz data/SHAPEIT/input.shapeit.22.gen.sample', 
                          'data/SHAPEIT/input.shapeit.22.hap.gz data/SHAPEIT/input.shapeit.22.hap.sample', 
                          'data/SHAPEIT/scaffold.haps.gz data/SHAPEIT/scaffold.haps.sample' ,'test')

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_runShapeit_woptions(clean_tmp):
    shapeit_folder=pytest.config.getoption("shapeit_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"
    work_dir= "data/outdir"

    command="perl {0}/standaloneJob.pl PyHive.VcfIntegration.run_Shapeit -language python3 \
    -shapeit_folder {1} -work_dir {2} -input_gen '{3}' -input_init '{4}' -input_scaffold '{5}' -outprefix {6} -inputthr 1.0 -thread 1 -window 0.1 \
    -states 400 -statesrandom 200 -burn 0 -run 12 -prune 4 -main 20 -inputfrom 20000000 -inputto 20100000 \
    -verbose True".format(hive_scripts, shapeit_folder, work_dir, 'data/SHAPEIT/input.shapeit.22.gen.gz data/SHAPEIT/input.shapeit.22.gen.sample',
                          'data/SHAPEIT/input.shapeit.22.hap.gz data/SHAPEIT/input.shapeit.22.hap.sample',
                          'data/SHAPEIT/scaffold.haps.gz data/SHAPEIT/scaffold.haps.sample' , 'test')

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)
