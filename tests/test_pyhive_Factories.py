import os
import subprocess
import glob
import pytest

# test_pyhive_Factories.py

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob("data/outdir/*")
    for f in files:
        os.remove(f)

def test_chrfactory():
    fa_ix = pytest.config.getoption("faix")
    hive_scripts = pytest.config.getoption("hive_lib")+"/scripts/"

    command = "perl {0}/standaloneJob.pl PyHive.Factories.ChrFactory " \
              "-language python3 -faix {1}".format(hive_scripts, fa_ix)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_BeagleChunkFactory():
    vcf_f = pytest.config.getoption("vcf_gts")
    beaglechunks_folder = pytest.config.getoption("makeBGLCHUNKS_folder")
    hive_scripts = pytest.config.getoption("hive_lib")+"/scripts/"
    work_dir = "data/outdir/"

    command = "perl {0}/standaloneJob.pl PyHive.Factories.BeagleChunkFactory " \
              "-language python3 -filepath {1} -makeBGLCHUNKS_folder {2} " \
              "-work_dir {3} -window {4} -overlap {5}".format(hive_scripts,
                                                              vcf_f, beaglechunks_folder,
                                                              work_dir, 100, 2)
    print(command)
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_bedtools_make_windows():
    bedtools_folder = pytest.config.getoption("bedtools_folder")
    hive_scripts = pytest.config.getoption("hive_lib")+"/scripts/"
    genome_f = 'data/chr1.genome'
    log_dir = 'data/outdir/'

    command = "perl {0}/standaloneJob.pl PyHive.Factories.CoordFactory -language python3 " \
              "-bedtools_folder {1} -genome_file {2} -window {3} -offsest {4} -log_dir {5}".\
        format(hive_scripts, bedtools_folder, genome_f, 100000000, 200000, log_dir)
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_bedtools_make_windows_w_chrom():
    bedtools_folder = pytest.config.getoption("bedtools_folder")
    hive_scripts = pytest.config.getoption("hive_lib")+"/scripts/"
    genome_f = 'data/chr1.genome'
    log_dir = 'data/outdir/'

    command = "perl {0}/standaloneJob.pl PyHive.Factories.CoordFactory " \
              "-language python3 -bedtools_folder {1} -genome_file {2} -window {3}\
    -offsest {4} -chrom chr1 -log_dir {5}".format(hive_scripts,
                                                  bedtools_folder,
                                                  genome_f,
                                                  100000000,
                                                  200000,
                                                  log_dir)
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_bedtools_make_windows_w_subtract(clean_tmp):
    '''
    Testing PyHive runnable with a BED file to subtract
    '''
    bedtools_folder = pytest.config.getoption("bedtools_folder")
    hive_scripts = pytest.config.getoption("hive_lib")+"/scripts/"
    genome_f = 'data/chr1.genome'
    subtract = 'data/subtract.bed'
    log_dir = 'data/outdir/'

    command = "perl {0}/standaloneJob.pl PyHive.Factories.CoordFactory -language python3 " \
              "-bedtools_folder {1} -genome_file {2} -window {3}\
    -offsest {4} -subtract {5} -log_dir {6}".format(hive_scripts, bedtools_folder,
                                                    genome_f, 100000000, 200000,
                                                    subtract, log_dir)
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)
