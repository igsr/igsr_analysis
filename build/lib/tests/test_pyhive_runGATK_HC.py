import os
import pytest
import subprocess
import glob

# test_pyhive_runGATK_HC.py

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/outdir/*')
    for f in files:
        os.remove(f)

def test_runGATK_HC(clean_tmp):

    shapeit_folder=pytest.config.getoption("shapeit_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"
    bam_file = pytest.config.getoption("--bam")
    reference = pytest.config.getoption("--reference")
    gatk_folder= pytest.config.getoption("--gatk_folder")
    bgzip_folder = pytest.config.getoption("--bgzip_folder")

    work_dir= "data/outdir/"

    command="perl {0}/standaloneJob.pl PyHive.VariantCalling.GATK_HC -language python3 \
    -outprefix {1} -work_dir {2} -chunk {3} -bamlist {4} -reference {5} \
    -gatk_folder {6} -bgzip_folder {7} -verbose True".format(hive_scripts, 'out', work_dir, 
                                                             "\"['chr1','10000','30000']\"", bam_file, 
                                                             reference, gatk_folder, bgzip_folder)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_runGATK_HC_with_alleles(clean_tmp):

    shapeit_folder=pytest.config.getoption("shapeit_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"
    bam_file = pytest.config.getoption("--bam")
    reference = pytest.config.getoption("--reference")
    gatk_folder= pytest.config.getoption("--gatk_folder")
    bgzip_folder = pytest.config.getoption("--bgzip_folder")
    alleles = pytest.config.getoption("--vcf")

    work_dir= "data/outdir/"

    command="perl {0}/standaloneJob.pl PyHive.VariantCalling.GATK_HC -language python3 \
    -outprefix {1} -work_dir {2} -chunk {3} -bamlist {4} -reference {5} \
    -gatk_folder {6} -bgzip_folder {7} -alleles {8} -verbose True".format(hive_scripts, 'out', work_dir,
                                                                          "\"['chr1','10000','30000']\"", bam_file,
                                                                          reference, gatk_folder, bgzip_folder, alleles)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_runGATK_HC_wlogile(clean_tmp):

    shapeit_folder=pytest.config.getoption("shapeit_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"
    bam_file = pytest.config.getoption("--bam")
    reference = pytest.config.getoption("--reference")
    gatk_folder= pytest.config.getoption("--gatk_folder")
    bgzip_folder = pytest.config.getoption("--bgzip_folder")

    work_dir= "data/outdir/"

    command="perl {0}/standaloneJob.pl PyHive.VariantCalling.GATK_HC -language python3 \
    -outprefix {1} -work_dir {2} -chunk {3} -bamlist {4} -reference {5} \
    -gatk_folder {6} -bgzip_folder {7} -log_file {8} -verbose True".format(hive_scripts, 'out', work_dir,
                                                                           "\"['chr1','10000','30000']\"", bam_file,
                                                                           reference, gatk_folder, bgzip_folder,
                                                                           "data/outdir/test_hc")
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_runGATK_HC_multithread(clean_tmp):

    shapeit_folder=pytest.config.getoption("shapeit_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"
    bam_file = pytest.config.getoption("--bam")
    reference = pytest.config.getoption("--reference")
    gatk_folder= pytest.config.getoption("--gatk_folder")
    bgzip_folder = pytest.config.getoption("--bgzip_folder")

    work_dir= "data/outdir/"

    command="perl {0}/standaloneJob.pl PyHive.VariantCalling.GATK_HC -language python3 \
    -outprefix {1} -work_dir {2} -chunk {3} -bamlist {4} -reference {5} \
    -gatk_folder {6} -bgzip_folder {7} -log_file {8} -threads 2 -verbose True".format(hive_scripts, 'out', work_dir,
                                                                                      "\"['chr1','10000','30000']\"", bam_file,
                                                                                      reference, gatk_folder, bgzip_folder,
                                                                                      "data/outdir/test_hc")
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)


