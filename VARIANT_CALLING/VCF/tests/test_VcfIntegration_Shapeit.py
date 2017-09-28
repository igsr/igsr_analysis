import os
import pytest
import glob
from VCFIntegration.Shapeit import Shapeit

# test_VcfIntegration_Shapeit.py

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/SHAPEIT/outdir/*')
    for f in files:
        os.remove(f)

def test_run_shapeit():

    shapeit_o=Shapeit(shapeit_folder = pytest.config.getoption("--shapeit_folder"))
    
    shapeit_o.run_shapeit(input_gen= 'data/SHAPEIT/input.shapeit.22.gen.gz data/SHAPEIT/input.shapeit.22.gen.sample',
                          input_init= 'data/SHAPEIT/input.shapeit.22.hap.gz data/SHAPEIT/input.shapeit.22.hap.sample',
                          input_scaffold= 'data/SHAPEIT/scaffold.haps.gz data/SHAPEIT/scaffold.haps.sample',
                          output_prefix='data/SHAPEIT/outdir/output.shapeit.22', verbose=True)
    assert os.path.exists('data/SHAPEIT/outdir/output.shapeit.22.haps.gz')

def test_run_shapeit_w_options(clean_tmp):

    shapeit_o=Shapeit(shapeit_folder = pytest.config.getoption("--shapeit_folder"))

    options={
        'input-from' : 20000000, 
        'input-to': 20100000,
        'input-thr' : 1.0,
        'thread' : 1,
        'window' : 0.1,
        'states' : 400,
        'states-random' : 200}

    shapeit_o.run_shapeit(input_gen= 'data/SHAPEIT/input.shapeit.22.gen.gz data/SHAPEIT/input.shapeit.22.gen.sample',
                          input_init= 'data/SHAPEIT/input.shapeit.22.hap.gz data/SHAPEIT/input.shapeit.22.hap.sample',
                          input_scaffold= 'data/SHAPEIT/scaffold.haps.gz data/SHAPEIT/scaffold.haps.sample',
                          output_prefix='data/SHAPEIT/outdir/output.shapeit.22.20000000.20100000',
                          verbose=True,
                          **options)
    assert os.path.exists('data/SHAPEIT/outdir/output.shapeit.22.20000000.20100000.haps.gz')
                          
