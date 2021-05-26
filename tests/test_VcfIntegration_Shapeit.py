import os
import glob
import pytest
import pdb
import shutil

from VCF.VCFIntegration import Shapeit

# test_VcfIntegration_Shapeit.py

@pytest.fixture
def shapeit_input(datadir):
    
    params = []
    params.append("{0}/SHAPEIT/input.shapeit.22.gen.gz {0}/SHAPEIT/input.shapeit.22.gen.sample".format(datadir))
    params.append("{0}/SHAPEIT/input.shapeit.22.hap.gz {0}/SHAPEIT/input.shapeit.22.hap.sample".format(datadir))
    params.append("{0}/SHAPEIT/scaffold.22.phased.haps {0}/SHAPEIT/scaffold.22.phased.sample".format(datadir))

    return params


def test_run_shapeit(shapeit_folder, shapeit_input, datadir):

    shapeit_o = Shapeit(shapeit_folder=shapeit_folder)

    shapeit_o.run_shapeit(input_gen=shapeit_input[0],
                          input_init=shapeit_input[1],
                          input_scaffold=shapeit_input[2],
                          output_prefix="{0}/outdir/output.shapeit.22".format(datadir), verbose=True)
    assert os.path.exists("{0}/outdir/output.shapeit.22.haps.gz".format(datadir))

def test_run_shapeit_w_options(shapeit_folder, shapeit_input, datadir):

    options = {
        'input-from': 20000000,
        'input-to': 20100000,
        'input-thr': 1.0,
        'thread': 1,
        'window': 0.1,
        'states': 400,
        'states-random': 200}

    shapeit_o = Shapeit(shapeit_folder=shapeit_folder)

    shapeit_o.run_shapeit(input_gen=shapeit_input[0],
                          input_init=shapeit_input[1],
                          input_scaffold=shapeit_input[2],
                          output_prefix="{0}/outdir/output.shapeit.22.20000000.20100000".format(datadir),
                          verbose=True,
                          **options)

    assert os.path.exists("{0}/outdir/output.shapeit.22.20000000.20100000.haps.gz".format(datadir))

def test_ligate_shapeitchunks(datadir):
    """
    Test method to run ligateHAPLOTYPES. This test throw an error because
     scaffolded_samples, chunk_str are fictitious
    """

    ligateHAPLOTYPES_folder = None # folder containing the makeBGLCHUNKS binary
    if shutil.which('ligateHAPLOTYPES') is None:
        raise Exception("'ligateHAPLOTYPES' needs to by in $PATH")
    else:
        ligateHAPLOTYPES_folder = os.path.dirname(shutil.which('ligateHAPLOTYPES'))

    shapeit_o = Shapeit(ligateHAPLOTYPES_folder=ligateHAPLOTYPES_folder)

    with pytest.raises(Exception):
        shapeit_o.ligate_shapeitchunks(vcf_f="{0}/test1.vcf.gz".format(datadir), scaffolded_samples='test.samples',
                                       chunk_str='s2.chunk1.hap.gz s2.chunk1.hap.gz s2.chunk1.hap.gz',
                                       output_prefix="{0}/outdir/test".format(datadir), verbose=True)

def test_run_shapeit_convert2vcf(shapeit_folder, datadir, clean_tmp):

    shapeit_o = Shapeit(shapeit_folder=shapeit_folder)

    out_vcf = shapeit_o.convert2vcf(input_prefix="{0}/outdir/output.shapeit.22.haps".format(datadir),
                                    output_prefix="{0}/outdir/output.shapeit.22.phased".format(datadir),
                                    compress=True,
                                    verbose=True,
                                    logfile="{0}/SHAPEIT.output.shapeit.22.phased.log".format(datadir))
    
    assert os.path.exists(out_vcf)
