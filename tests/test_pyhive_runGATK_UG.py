import os
import pytest
import subprocess
import glob
import pdb

# test_pyhive_runGATK_UG.py

def test_runGATK_UG(hive_dir, bgzip_folder, gatk_folder, datadir, clean_tmp):

    bam_file = "{0}/exampleBAM.bam".format(datadir)
    reference = "{0}/exampleFASTA.fasta".format(datadir)
    work_dir = "{0}/outdir".format(datadir)

    command="perl {0}/scripts/standaloneJob.pl PyHive.VariantCalling.GATK_UG -language python3 \
    -outprefix {1} -work_dir {2} -chunk {3} -bamlist {4} -reference {5} \
    -gatk_folder {6} -bgzip_folder {7} -glm SNP -output_mode EMIT_ALL_SITES -verbose True".format(hive_dir, 'out', work_dir, 
                                                                                       "\"['chr1','10000','30000']\"", bam_file, 
                                                                                       reference, gatk_folder, bgzip_folder)
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_runGATK_UG_with_alleles(hive_dir, bgzip_folder, gatk_folder, datadir, clean_tmp):

    bam_file = "{0}/exampleBAM.bam".format(datadir)
    reference = "{0}/exampleFASTA.fasta".format(datadir)
    alleles = "{0}/test.vcf.gz".format(datadir)
    work_dir = "{0}/outdir".format(datadir)

    command="perl {0}/scripts/standaloneJob.pl PyHive.VariantCalling.GATK_UG -language python3 \
    -outprefix {1} -work_dir {2} -chunk {3} -bamlist {4} -reference {5} \
    -gatk_folder {6} -bgzip_folder {7} -alleles {8} -glm SNP -output_mode EMIT_ALL_SITES -verbose True".format(hive_dir, 'out', work_dir,
                                                                                                     "\"['chr1','10000','30000']\"", bam_file,
                                                                                                     reference, gatk_folder, bgzip_folder, alleles)
    
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_runGATK_UG_wlogile(hive_dir, bgzip_folder, gatk_folder, datadir, clean_tmp):

    bam_file = "{0}/exampleBAM.bam".format(datadir)
    reference = "{0}/exampleFASTA.fasta".format(datadir)
    work_dir = "{0}/outdir".format(datadir)

    command="perl {0}/scripts/standaloneJob.pl PyHive.VariantCalling.GATK_UG -language python3 \
    -outprefix {1} -work_dir {2} -chunk {3} -bamlist {4} -reference {5} \
    -gatk_folder {6} -bgzip_folder {7} -glm SNP -output_mode EMIT_ALL_SITES -log_file {8} -verbose True".format(hive_dir, 'out', work_dir,
                                                                                                      "\"['chr1','10000','30000']\"", bam_file,
                                                                                                      reference, gatk_folder, bgzip_folder,
                                                                                                      "data/outdir/test")
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_runGATK_UG_throws_exception(hive_dir, bgzip_folder, gatk_folder, datadir, clean_tmp):
    """
    Test eHive runnable for GATK UnifiedGenotyper and throws an exception
    """

    bam_file = "{0}/exampleBAM.bam".format(datadir)
    reference = "{0}/exampleFASTA.fasta".format(datadir)
    work_dir = "{0}/outdir".format(datadir)
    glm = 'fake'
    output_mode = 'EMIT_ALL_SITES'

    command="perl {0}/scripts/standaloneJob.pl PyHive.VariantCalling.GATK_UG -language python3 \
    -outprefix {1} -work_dir {2} -chunk {3} -bamlist {4} -reference {5} \
    -gatk_folder {6} -bgzip_folder {7} -glm {8} -output_mode {9} -verbose True".format(hive_dir, 'out', work_dir,
                                                                                       "\"['chr1','10000','30000']\"", bam_file,
                                                                                       reference, gatk_folder, bgzip_folder,
                                                                                       glm, output_mode)
    with pytest.raises(Exception):
        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            raise Exception(exc.output)
