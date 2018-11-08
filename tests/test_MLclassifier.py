import os
import warnings
import pytest
import glob

from VCF.VCFfilter.MLclassifier import MLclassifier

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/outdir/*')
    for f in files:
        os.remove(f)

def test_train_snps():
    '''
    Train the model for SNPs
    '''

    ML_obj=MLclassifier(bcftools_folder = pytest.config.getoption('bcftools_folder'))
    outfile=ML_obj.train(outprefix="data/outdir/fitted_logreg_snps", 
                         tp_annotations=pytest.config.getoption('--tp_annotations_snps'),
                         fp_annotations=pytest.config.getoption('--fp_annotations_snps'))
   
    assert os.path.isfile(outfile) is True

def test_train_indels(clean_tmp):
    '''
    Train the model for INDELs
    '''

    ML_obj=MLclassifier(bcftools_folder = pytest.config.getoption('bcftools_folder'))
    outfile=ML_obj.train(outprefix="data/outdir/fitted_logreg_indels",
                         tp_annotations=pytest.config.getoption('--tp_annotations_indels'),
                         fp_annotations=pytest.config.getoption('--fp_annotations_indels'))

    assert os.path.isfile(outfile) is True
