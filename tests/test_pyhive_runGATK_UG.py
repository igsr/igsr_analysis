import os
import pytest
import subprocess
import glob

# test_pyhive_runGATK_UG.py

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/outdir/*')
    for f in files:
        os.remove(f)

<<<<<<< HEAD

=======
"""
>>>>>>> origin/newphasing
def test_runGATK_UG(clean_tmp):

    shapeit_folder=pytest.config.getoption("shapeit_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"
    bam_file = pytest.config.getoption("--bam")
    reference = pytest.config.getoption("--reference")
    gatk_folder= pytest.config.getoption("--gatk_folder")
    bgzip_folder = pytest.config.getoption("--bgzip_folder")
    glm = pytest.config.getoption("--glm")
    output_mode = pytest.config.getoption("--output_mode")

    work_dir= "data/outdir/"

    command="perl {0}/standaloneJob.pl PyHive.VariantCalling.GATK_UG -language python3 \
    -outprefix {1} -work_dir {2} -chunk {3} -bamlist {4} -reference {5} \
    -gatk_folder {6} -bgzip_folder {7} -glm {8} -output_mode {9} -verbose True".format(hive_scripts, 'out', work_dir, 
                                                                                       "\"['chr1','10000','30000']\"", bam_file, 
                                                                                       reference, gatk_folder, bgzip_folder,
                                                                                       glm, output_mode)
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)
"""

def test_runGATK_UG_with_alleles(clean_tmp):

    shapeit_folder=pytest.config.getoption("shapeit_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"
    bam_file = pytest.config.getoption("--bam")
    reference = pytest.config.getoption("--reference")
    gatk_folder= pytest.config.getoption("--gatk_folder")
    bgzip_folder = pytest.config.getoption("--bgzip_folder")
    alleles = pytest.config.getoption("--vcf")
    glm = pytest.config.getoption("--glm")
    output_mode = pytest.config.getoption("--output_mode")

    work_dir= "data/outdir/"

    command="perl {0}/standaloneJob.pl PyHive.VariantCalling.GATK_UG -language python3 \
    -outprefix {1} -work_dir {2} -chunk {3} -bamlist {4} -reference {5} \
    -gatk_folder {6} -bgzip_folder {7} -alleles {8} -glm {9} -output_mode {10} -verbose True".format(hive_scripts, 'out', work_dir,
                                                                                                     "\"['chr1','10000','30000']\"", bam_file,
                                                                                                     reference, gatk_folder, bgzip_folder, alleles,
                                                                                                     glm, output_mode)

<<<<<<< HEAD
=======
    print(command)
>>>>>>> origin/newphasing

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

"""
def test_runGATK_UG_wlogile(clean_tmp):

    shapeit_folder=pytest.config.getoption("shapeit_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"
    bam_file = pytest.config.getoption("--bam")
    reference = pytest.config.getoption("--reference")
    gatk_folder= pytest.config.getoption("--gatk_folder")
    bgzip_folder = pytest.config.getoption("--bgzip_folder")
    glm = pytest.config.getoption("--glm")
    output_mode = pytest.config.getoption("--output_mode")

    work_dir= "data/outdir/"

    command="perl {0}/standaloneJob.pl PyHive.VariantCalling.GATK_UG -language python3 \
    -outprefix {1} -work_dir {2} -chunk {3} -bamlist {4} -reference {5} \
    -gatk_folder {6} -bgzip_folder {7} -glm {8} -output_mode {9} -log_file {10} -verbose True".format(hive_scripts, 'out', work_dir,
                                                                                                      "\"['chr1','10000','30000']\"", bam_file,
                                                                                                      reference, gatk_folder, bgzip_folder,
                                                                                                      glm, output_mode, "data/outdir/test")
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_runGATK_UG_throws_exception(clean_tmp):
    '''
    Test eHive runnable for GATK UnifiedGenotyper and throws an exception
    '''

    shapeit_folder=pytest.config.getoption("shapeit_folder")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/"
    bam_file = pytest.config.getoption("--bam")
    reference = pytest.config.getoption("--reference")
    gatk_folder= pytest.config.getoption("--gatk_folder")
    bgzip_folder = pytest.config.getoption("--bgzip_folder")
    glm = 'fake'
    output_mode = pytest.config.getoption("--output_mode")

    work_dir= "data/outdir/"

    command="perl {0}/standaloneJob.pl PyHive.VariantCalling.GATK_UG -language python3 \
    -outprefix {1} -work_dir {2} -chunk {3} -bamlist {4} -reference {5} \
    -gatk_folder {6} -bgzip_folder {7} -glm {8} -output_mode {9} -verbose True".format(hive_scripts, 'out', work_dir,
                                                                                       "\"['chr1','10000','30000']\"", bam_file,
                                                                                       reference, gatk_folder, bgzip_folder,
                                                                                       glm, output_mode)
    with pytest.raises(Exception):
        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            raise Exception(exc.output)

<<<<<<< HEAD
=======
"""
>>>>>>> origin/newphasing
