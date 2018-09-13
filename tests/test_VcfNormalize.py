import os
import pytest
import glob
import pdb

from VCF.VcfNormalize import VcfNormalize

# test_VcfNormalize.py

@pytest.fixture
def vcf_object():
    '''Returns an  object'''

    vcf_file = pytest.config.getoption("--vcf")
    bcftools_folder = pytest.config.getoption("--bcftools_folder")
    bgzip_folder = pytest.config.getoption("--bgzip_folder")
    vcflib_folder = pytest.config.getoption("--vcflib_folder")
    vt_folder = pytest.config.getoption("--vt_folder")
    gatk_folder = pytest.config.getoption("--gatk_folder")

    vcf_object=VcfNormalize(vcf=vcf_file,bgzip_folder=bgzip_folder,
                            vcflib_folder=vcflib_folder, vt_folder=vt_folder,
                            bcftools_folder=bcftools_folder,gatk_folder=gatk_folder)
    return vcf_object

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/outdir/*')
    for f in files:
        os.remove(f)
"""
def test_run_vtnormalize(vcf_object):
    outfile=vcf_object.run_vtnormalize(outprefix='data/outdir/test', reference=pytest.config.getoption("--reference"),verbose=True, n=True)
    
    assert os.path.exists(outfile)

def test_run_vtnormalize_compress(vcf_object):
    '''
    Test vt normalize with compress=True option
    '''
    outfile=vcf_object.run_vtnormalize(outprefix='data/outdir/test', reference=pytest.config.getoption("--reference"),verbose=True, n=True,compress=True)

    assert os.path.exists(outfile)


def test_run_vcfallelicprimitives(vcf_object):
    outfile=vcf_object.run_vcfallelicprimitives(outprefix='data/outdir/test',verbose=True)

    assert os.path.exists(outfile)

def test_run_vcfallelicprimitives_downstream_pipe(vcf_object):
    '''
    Test function to run vcfallelicprimitives and piping to other programs
    '''
    outfile=vcf_object.run_vcfallelicprimitives(outprefix='data/outdir/test1', 
                                                downstream_pipe='~/bin/vt/vt sort -',
                                                verbose=True, compress=True)
    assert os.path.exists(outfile)
"""
def test_run_GATK_VariantsToAllelicPrimitives(vcf_object, clean_tmp):
    '''
    Test function to run GATK VariantsToAllelicPrimitives
    '''

    outfile=vcf_object.run_gatk_VariantsToAllelicPrimitives(outprefix='data/outdir/test2',
                                                            reference=pytest.config.getoption("--reference"),
                                                            compress=True, verbose=True)
    assert os.path.exists(outfile)
   
