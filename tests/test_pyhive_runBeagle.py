import os
import subprocess
import glob
import pytest
import pdb

from PyHive.VcfIntegration import *

# test_pyhive_runBeagle.py

def test_runBeagle(beagle_bins, hive_dir, datadir, clean_tmp):

    vcf_f = "{0}/test_chr20.vcf.gz".format(datadir)
    beagle_folder, beagle_jar = beagle_bins

    work_dir = "{0}/outdir".format(datadir)

    command = "perl {0}/scripts/standaloneJob.pl PyHive.VcfIntegration.run_Beagle -language python3 " \
              "-vcf_file {1} -beagle_folder {2} -beagle_jar {3} -work_dir {4} -outprefix {5} " \
              "-window {6} -overlap {7} -niterations {8} -nthreads {9} -region_chunk {10} " \
              "-verbose True".format(hive_dir, vcf_f, beagle_folder, beagle_jar, work_dir,
                                     'test', 12000, 2000, 15, 1, 'chr20:1000000-1000150')
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_runBeagle_noparams(beagle_bins, hive_dir, datadir, clean_tmp):

    vcf_f = "{0}/test_chr20.vcf.gz".format(datadir)
    beagle_folder, beagle_jar = beagle_bins
    
    work_dir = "{0}/outdir".format(datadir)

    command = "perl {0}/scripts/standaloneJob.pl PyHive.VcfIntegration.run_Beagle -language python3 " \
              "-vcf_file {1} -beagle_folder {2} -beagle_jar {3} -work_dir {4} -outprefix {5} " \
              "-region_chunk {6} -verbose True".format(hive_dir, vcf_f, beagle_folder,
                                                       beagle_jar,
                                                       work_dir,
                                                       'test',
                                                       'chr20:1000000-1000150')

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_runBeagle_noparams_correct(beagle_bins, hive_dir, datadir, clean_tmp):
    
    vcf_f = "{0}/test_chr20.vcf.gz".format(datadir)
    beagle_folder, beagle_jar = beagle_bins
    
    work_dir = "{0}/outdir".format(datadir)

    command = "perl {0}/scripts/standaloneJob.pl PyHive.VcfIntegration.run_Beagle -language python3 " \
              "-vcf_file {1} -beagle_folder {2} -beagle_jar {3} -work_dir {4} -outprefix {5} " \
              "-region_chunk {6} -correct True -verbose True".format(hive_dir,
                                                                     vcf_f, beagle_folder,
                                                                     beagle_jar, work_dir,
                                                                     'test_correct',
                                                                     'chr20:1000000-1000150')
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)
