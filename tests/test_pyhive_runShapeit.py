import os
import subprocess
import glob
import pytest
import pdb

from PyHive.VcfIntegration import *

# test_pyhive_runShapeit.py

def test_runShapeit(hive_dir, shapeit_folder, datadir, clean_tmp):

    work_dir = "{0}/outdir".format(datadir)

    shapeit_input1 = "{0}/SHAPEIT/input.shapeit.22.gen.gz {0}/SHAPEIT/input.shapeit.22.gen.sample".format(datadir)
    shapeit_input2 = "{0}/SHAPEIT/input.shapeit.22.hap.gz {0}/SHAPEIT/input.shapeit.22.hap.sample".format(datadir)
    shapeit_input3 = "{0}/SHAPEIT/scaffold".format(datadir)

    command = "perl {0}/scripts/standaloneJob.pl PyHive.VcfIntegration.run_Shapeit -language python3 \
    -shapeit_folder {1} -work_dir {2} -input_gen '{3}' -input_init '{4}' -input_scaffold_prefix" \
              " \"['{5}']\" -chr chr22 -outprefix {6} -verbose True".\
        format(hive_dir, shapeit_folder, work_dir, shapeit_input1,
               shapeit_input2, shapeit_input3, 'test')

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_runShapeit_woptions(hive_dir, shapeit_folder, datadir, clean_tmp):

    work_dir = "{0}/outdir".format(datadir)

    shapeit_input1 = "{0}/SHAPEIT/input.shapeit.22.gen.gz {0}/SHAPEIT/input.shapeit.22.gen.sample".format(datadir)
    shapeit_input2 = "{0}/SHAPEIT/input.shapeit.22.hap.gz {0}/SHAPEIT/input.shapeit.22.hap.sample".format(datadir)
    shapeit_input3 = "{0}/SHAPEIT/scaffold".format(datadir)

    command = "perl {0}/scripts/standaloneJob.pl PyHive.VcfIntegration.run_Shapeit -language python3 \
    -shapeit_folder {1} -work_dir {2} -input_gen '{3}' -input_init '{4}' -input_scaffold_prefix" \
              " \"['{5}']\" -outprefix {6} -chr chr22 -inputthr 1.0 -thread 1 -window 0.1 \
    -states 400 -statesrandom 200 -burn 0 -run 12 -prune 4 -main 20" \
              " -inputfrom 20000000 -inputto 20100000 \
    -verbose True".format(hive_dir, shapeit_folder, work_dir, shapeit_input1, shapeit_input2,
                           shapeit_input3, 'test')

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)

def test_run_Shapeit_convert2vcf(hive_dir, shapeit_folder, datadir, clean_tmp):

    work_dir = "{0}/outdir".format(datadir)

    command = "perl {0}/scripts/standaloneJob.pl PyHive.VcfIntegration." \
              "run_Shapeit_convert2vcf -language python3 \
    -shapeit_folder {1} -work_dir {2} -hap_gz {3} -hap_sample {4} -outprefix {5} \
    -compress True -verbose True".format(hive_dir,
                                         shapeit_folder,
                                         work_dir,
                                         "{0}/SHAPEIT/input.shapeit.22.hap.gz".format(datadir),
                                         "{0}/SHAPEIT/input.shapeit.22.hap.sample".format(datadir),
                                         'test.phased')

    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)
