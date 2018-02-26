import os
import pytest
import subprocess
from BamQC import BamQC

# test_pyhive_BamQC.py

def test_RunChkIndelRg():
    fa_ix= pytest.config.getoption("bam_f")
    hive_scripts= pytest.config.getoption("hive_lib")+"/scripts/" 

    command="perl {0}/standaloneJob.pl PyHive.Factories.ChrFactory -language python3 -faix {1}".format(hive_scripts, fa_ix)
    
    try:
        subprocess.check_output(command, shell=True)
        assert True
    except subprocess.CalledProcessError as exc:
        assert False
        raise Exception(exc.output)
    
