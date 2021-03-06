import os
import subprocess
import glob
import pytest

# test_pyhive_runBCFTools_VC.py

def test_runBCFTools_VC(bcftools_folder, hive_dir, datadir, clean_tmp):
    """
    Test function to run BCFTools mpileup|call on a BAM file
    """

    bam_file = "{0}/exampleBAM.bam".format(datadir)
    reference = "{0}/exampleFASTA.fasta".format(datadir)

    work_dir = "{0}/outdir/".format(datadir)
    annots = "\"['DP','SP','AD']\""

    command = "perl {0}/scripts/standaloneJob.pl PyHive.VariantCalling.BCFTools_caller -language python3 \
    -outprefix {1} -work_dir {2} -chunk {3} -bam {4} -reference {5} \
    -bcftools_folder {6} -annots {7} -verbose True".format(hive_dir, 'out', work_dir,
                                                           "\"['chr1','10000','30000']\"", bam_file,
                                                           reference, bcftools_folder, annots)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_runBCFTools_VC_woptions(bcftools_folder, hive_dir, datadir, clean_tmp):
    """
    Test function to run BCFTools mpileup|call on a BAM file
    using some options and arguments
    """

    bam_file = "{0}/exampleBAM.bam".format(datadir)
    reference = "{0}/exampleFASTA.fasta".format(datadir)

    work_dir = "{0}/outdir/".format(datadir)
    annots = "\"['DP','SP','AD']\""

    command = "perl {0}/scripts/standaloneJob.pl PyHive.VariantCalling.BCFTools_caller -language python3 \
    -outprefix {1} -work_dir {2} -chunk {3} -bam {4} -reference {5} \
    -bcftools_folder {6} -annots {7} -E 1 -p 1 -m_pileup 3 -m_call 1 -v 1 " \
              "-F 0.05 -C 25 -verbose True".format(hive_dir, 'out', work_dir,
                                                   "\"['chr1','10000','30000']\"", bam_file,
                                                   reference, bcftools_folder, annots)

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)
