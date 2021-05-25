import os
import subprocess
import glob
import pytest
import pdb

# test_pyhive_runGATK_HC.py

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/outdir/*')
    for f in files:
        os.remove(f)

def test_runGATK_HC(hive_dir, bgzip_folder, gatk_folder, datadir, clean_tmp):

    bam_file = "{0}/exampleBAM.bam".format(datadir)
    reference = "{0}/exampleFASTA.fasta".format(datadir)
    work_dir = "{0}/outdir".format(datadir)

    command = "perl {0}/scripts/standaloneJob.pl PyHive.VariantCalling.GATK_HC -language python3 \
    -outprefix {1} -work_dir {2} -chunk {3} -bamlist {4} -reference {5} \
    -gatk_folder {6} -bgzip_folder {7} -verbose True".\
        format(hive_dir, 'out', work_dir,
               "\"['chr1','10000','30000']\"",
               bam_file, reference, gatk_folder, bgzip_folder)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_runGATK_HC_with_alleles(hive_dir, bgzip_folder, gatk_folder, datadir, clean_tmp):

    bam_file = "{0}/exampleBAM.bam".format(datadir)
    reference = "{0}/exampleFASTA.fasta".format(datadir)

    alleles = "{0}/test.vcf.gz".format(datadir)

    work_dir = "{0}/outdir".format(datadir)

    command = "perl {0}/scripts/standaloneJob.pl PyHive.VariantCalling.GATK_HC -language python3 \
    -outprefix {1} -work_dir {2} -chunk {3} -bamlist {4} -reference {5} \
    -gatk_folder {6} -bgzip_folder {7} -alleles {8} " \
              "-verbose True".format(hive_dir, 'out', work_dir,
                                     "\"['chr1','10000','30000']\"", bam_file,
                                     reference, gatk_folder, bgzip_folder, alleles)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_runGATK_HC_wlogile(hive_dir, bgzip_folder, gatk_folder, datadir, clean_tmp):

    bam_file = "{0}/exampleBAM.bam".format(datadir)
    reference = "{0}/exampleFASTA.fasta".format(datadir)

    work_dir = "{0}/outdir".format(datadir)

    command = "perl {0}/scripts/standaloneJob.pl PyHive.VariantCalling.GATK_HC " \
              "-language python3 -outprefix {1} -work_dir {2} -chunk {3} " \
              "-bamlist {4} -reference {5} -gatk_folder {6} -bgzip_folder {7} " \
              "-log_file {8} -verbose True".format(hive_dir, 'out', work_dir,
                                                   "\"['chr1','10000','30000']\"", bam_file,
                                                   reference, gatk_folder, bgzip_folder,
                                                   "data/outdir/test_hc")
    
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

