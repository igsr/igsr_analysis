import os
import subprocess
import glob
import pytest

from PyHive.VcfIntegration import *

# test_pyhive_runShapeit.py

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/outdir/*')
    for f in files:
        os.remove(f)

def test_runShapeit():
    shapeit_folder = pytest.config.getoption("shapeit_folder")
    hive_scripts = pytest.config.getoption("hive_lib")+"/scripts/"
    work_dir = "data/outdir"

    command = "perl {0}/standaloneJob.pl PyHive.VcfIntegration.run_Shapeit -language python3 \
    -shapeit_folder {1} -work_dir {2} -input_gen '{3}' -input_init '{4}' -input_scaffold_prefix" \
              " \"['{5}']\" -chr chr22 -outprefix {6} -verbose True".\
        format(hive_scripts, shapeit_folder, work_dir,
               'data/SHAPEIT/input.shapeit.22.gen.gz '
               'data/SHAPEIT/input.shapeit.22.gen.sample',
               'data/SHAPEIT/input.shapeit.22.hap.gz'
               'data/SHAPEIT/input.shapeit.22.hap.sample',
               'data/SHAPEIT/scaffold',
               'test')
    print(command)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_runShapeit_woptions():
    shapeit_folder = pytest.config.getoption("shapeit_folder")
    hive_scripts = pytest.config.getoption("hive_lib")+"/scripts/"
    work_dir = "data/outdir"

    command = "perl {0}/standaloneJob.pl PyHive.VcfIntegration.run_Shapeit -language python3 \
    -shapeit_folder {1} -work_dir {2} -input_gen '{3}' -input_init '{4}' -input_scaffold_prefix" \
              " \"['{5}']\" -outprefix {6} -chr chr22 -inputthr 1.0 -thread 1 -window 0.1 \
    -states 400 -statesrandom 200 -burn 0 -run 12 -prune 4 -main 20" \
              " -inputfrom 20000000 -inputto 20100000 \
    -verbose True".format(hive_scripts, shapeit_folder, work_dir,
                          'data/SHAPEIT/input.shapeit.22.gen.gz '
                          'data/SHAPEIT/input.shapeit.22.gen.sample',
                          'data/SHAPEIT/input.shapeit.22.hap.gz '
                          'data/SHAPEIT/input.shapeit.22.hap.sample',
                          'data/SHAPEIT/scaffold',
                          'test')

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_run_Shapeit_convert2vcf(clean_tmp):
    shapeit_folder = pytest.config.getoption("shapeit_folder")
    hive_scripts = pytest.config.getoption("hive_lib")+"/scripts/"
    work_dir = "data/outdir"

    command = "perl {0}/standaloneJob.pl PyHive.VcfIntegration." \
              "run_Shapeit_convert2vcf -language python3 \
    -shapeit_folder {1} -work_dir {2} -hap_gz {3} -hap_sample {4} -outprefix {5} \
    -compress True -verbose True".format(hive_scripts,
                                         shapeit_folder,
                                         work_dir,
                                         'data/SHAPEIT/input.shapeit.22.hap.gz',
                                         'data/SHAPEIT/input.shapeit.22.hap.sample',
                                         'test.phased')

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)
