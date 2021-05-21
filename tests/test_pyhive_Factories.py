import os
import subprocess
import glob
import pytest
import pdb
import shutil

# test_pyhive_Factories.py

def test_chrfactory(hive_dir, datadir):

    fa_ix = "{0}/canonical_chros.fa.fai".format(datadir)

    command = "perl {0}/scripts/standaloneJob.pl PyHive.Factories.ChrFactory " \
              "-language python3 -faix {1}".format(hive_dir, fa_ix)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_BeagleChunkFactory(hive_dir, clean_tmp, datadir):

    vcf_f = "{0}/GLs.HG00136.vcf.gz".format(datadir)
    beaglechunks_folder = None
    if shutil.which('makeBGLCHUNKS') is None:
        raise Exception("'makeBGLCHUNKS needs to by in $PATH")
    else:
        beaglechunks_folder = os.path.dirname(shutil.which('makeBGLCHUNKS'))

    work_dir = "{0}/outdir/".format(datadir)

    command = "perl {0}/scripts/standaloneJob.pl PyHive.Factories.BeagleChunkFactory " \
              "-language python3 -filepath {1} -makeBGLCHUNKS_folder {2} " \
              "-work_dir {3} -window {4} -overlap {5}".format(hive_dir,
                                                              vcf_f, beaglechunks_folder,
                                                              work_dir, 100, 2)
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_bedtools_make_windows(bedtools_folder, hive_dir, datadir):

    genome_f = "{0}/chr1.genome".format(datadir)
    log_dir = "{0}/outdir/".format(datadir)

    command = "perl {0}/scripts/standaloneJob.pl PyHive.Factories.CoordFactory -language python3 " \
              "-bedtools_folder {1} -genome_file {2} -window {3} -offsest {4} -log_dir {5}".\
        format(hive_dir, bedtools_folder, genome_f, 100000000, 200000, log_dir)
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_bedtools_make_windows_w_chrom(bedtools_folder, hive_dir, datadir, clean_tmp):

    genome_f = "{0}/chr1.genome".format(datadir)
    log_dir = "{0}/outdir/".format(datadir)

    command = "perl {0}/scripts/standaloneJob.pl PyHive.Factories.CoordFactory " \
              "-language python3 -bedtools_folder {1} -genome_file {2} -window {3}\
    -offsest {4} -chrom chr1 -log_dir {5}".format(hive_dir,
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

def test_bedtools_make_windows_w_subtract(bedtools_folder, hive_dir, datadir, clean_tmp):
    """
    Testing PyHive runnable with a BED file to subtract
    """

    genome_f = "{0}/chr1.genome".format(datadir)
    subtract = "{0}/subtract.bed".format(datadir)
    log_dir = "{0}/outdir/".format(datadir)

    command = "perl {0}/scripts/standaloneJob.pl PyHive.Factories.CoordFactory -language python3 " \
              "-bedtools_folder {1} -genome_file {2} -window {3} "\
              "-offsest {4} -subtract {5} -log_dir {6}".format(hive_dir, bedtools_folder, genome_f, 100000000, 200000,
                                                    subtract, log_dir)
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)
