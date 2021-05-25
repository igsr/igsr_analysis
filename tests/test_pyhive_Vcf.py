import os
import subprocess
import glob
import pytest
import pdb

#from VcfUtils import VcfUtils

def test_VcfReheader(datadir, bcftools_folder, hive_dir, clean_tmp):
    """
    Test PyHive.Vcf.VcfReheader after setting the 'work_dir' arg
    """
    vcf_file = "{0}/test.vcf.gz".format(datadir)
    work_dir = "{0}/outdir/".format(datadir)

    command = "perl {0}/scripts/standaloneJob.pl PyHive.Vcf.VcfReheader " \
              "-language python3 -filepath {1} -bcftools_folder {2}" \
              " -newheader data/newheader.txt -work_dir {3} -debug 1".\
        format(hive_dir, vcf_file, bcftools_folder, work_dir)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_VcfReheader_w_samplename(datadir, bcftools_folder, hive_dir, clean_tmp):
    """
    Test PyHive.Vcf.VcfReheader after setting the 'samplename' arg
    """
    vcf_file = "{0}/test.vcf.gz".format(datadir)
    work_dir = "{0}/outdir/".format(datadir)

    command = "perl {0}/scripts/standaloneJob.pl PyHive.Vcf.VcfReheader -language python3 " \
              "-filepath {1} -bcftools_folder {2} -newheader data/newheader.txt " \
              "-work_dir {3} -samplename testname -debug 1".format(hive_dir,
                                                                    vcf_file,
                                                                    bcftools_folder,
                                                                    work_dir)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_VcfReheader_w_samplenfile(datadir, bcftools_folder, hive_dir, clean_tmp):
    """
    Test PyHive.Vcf.VcfReheader after setting the 'samplename' arg
    """
    vcf_file = "{0}/test.vcf.gz".format(datadir)
    work_dir = "{0}/outdir/".format(datadir)

    command = "perl {0}/scripts/standaloneJob.pl PyHive.Vcf.VcfReheader -language python3 -filepath {1} " \
              "-bcftools_folder {2} -newheader data/newheader.txt -work_dir {3} " \
              "-samplefile data/samplefile.txt -debug 1".format(hive_dir,
                                                                vcf_file,
                                                                bcftools_folder,
                                                                work_dir)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_VcfReplaceChrNames(datadir, bgzip_folder, hive_dir, clean_tmp):
    """
    Test PyHive.Vcf.VcfReplaceChrNames
    """
    vcf_file = "{0}/test.vcf.gz".format(datadir)
    work_dir = "{0}/outdir/".format(datadir)

    command = "perl {0}/scripts/standaloneJob.pl PyHive.Vcf.VcfReplaceChrNames -language python3 " \
              "-filepath {1} -chr_types 'ensembl' -work_dir {2} -outprefix test.ensembl " \
              "-bgzip_folder {3}".format(hive_dir, vcf_file, work_dir, bgzip_folder)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_VcfconvertPL2GL(datadir, bcftools_folder, hive_dir, clean_tmp):
    '''
    Test PyHive.Vcf.convertPL2GL
    '''
    vcf_file = "{0}/test.gatk.vcf.gz".format(datadir)
    work_dir = "{0}/outdir/".format(datadir)
   
    command = "perl {0}/scripts/standaloneJob.pl PyHive.Vcf.convertPL2GL -language python3 " \
              "-filepath {1} -work_dir {2} -outprefix test.pl2gl " \
              "-bcftools_folder {3}".format(hive_dir, vcf_file, work_dir, bcftools_folder)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)
