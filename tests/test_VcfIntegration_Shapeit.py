import os
import pytest
import glob
from VCFIntegration import Shapeit

# test_VcfIntegration_Shapeit.py

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/outdir/*')
    for f in files:
        os.remove(f)

def test_run_shapeit():

    shapeit_o=Shapeit(shapeit_folder = pytest.config.getoption("--shapeit_folder"))
    
    shapeit_o.run_shapeit(input_gen= 'data/SHAPEIT/input.shapeit.22.gen.gz data/SHAPEIT/input.shapeit.22.gen.sample',
                          input_init= 'data/SHAPEIT/input.shapeit.22.hap.gz data/SHAPEIT/input.shapeit.22.hap.sample',
                          input_scaffold= 'data/SHAPEIT/scaffold.22.phased.haps data/SHAPEIT/scaffold.22.phased.sample',
                          output_prefix='data/outdir/output.shapeit.22', verbose=True)
    assert os.path.exists('data/outdir/output.shapeit.22.haps.gz')

def test_run_shapeit_w_options():

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
                          input_scaffold= 'data/SHAPEIT/scaffold.22.phased.haps data/SHAPEIT/scaffold.22.phased.sample',
                          output_prefix='data/outdir/output.shapeit.22.20000000.20100000',
                          verbose=True,
                          **options)
    assert os.path.exists('data/outdir/output.shapeit.22.20000000.20100000.haps.gz')

def test_ligate_shapeitchunks():
    '''
    Test method to run ligateHAPLOTYPES. This test throw an error because scaffolded_samples, chunk_str  are ficticious
    '''

    shapeit_o=Shapeit(ligateHAPLOTYPES_folder = pytest.config.getoption("--ligateHAPLOTYPES_folder"))

    with pytest.raises(Exception):
        shapeit_o.ligate_shapeitchunks(vcf_f='data/test1.vcf.gz',scaffolded_samples='test.samples',
                                       chunk_str='s2.chunk1.hap.gz s2.chunk1.hap.gz s2.chunk1.hap.gz',
                                       output_prefix='data/outdir/test',verbose=True)

def test_run_shapeit_convert2vcf(clean_tmp):

    shapeit_o=Shapeit(shapeit_folder = pytest.config.getoption("--shapeit_folder"))

    out_vcf=shapeit_o.convert2vcf(input_prefix='data/outdir/output.shapeit.22.haps', 
                                  output_prefix='data/outdir/output.shapeit.22.phased', 
                                  compress=True, verbose=True, logfile='data/SHAPEIT.output.shapeit.22.phased.log')
    
    assert os.path.exists(out_vcf)

