import os
import glob
import pytest

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

    ML_obj = MLclassifier(bcftools_folder=pytest.config.getoption('bcftools_folder'))
    outfile = ML_obj.train(outprefix="data/outdir/fitted_logreg_snps",
                           tp_annotations=pytest.config.getoption('--tp_annotations_snps'),
                           fp_annotations=pytest.config.getoption('--fp_annotations_snps'))

    assert os.path.isfile(outfile) is True

def test_train_snps_gz():
    '''
    Train the model for SNPs using 2 gzipped annotation files
    '''

    ML_obj = MLclassifier(bcftools_folder=pytest.config.getoption('bcftools_folder'))
    outfile = ML_obj.train(outprefix="data/outdir/fitted_logreg_snps",
                           tp_annotations=pytest.config.getoption('--tp_annotations_snps_gz'),
                           fp_annotations=pytest.config.getoption('--fp_annotations_snps_gz'))

def test_train_indels():
    '''
    Train the model for INDELs
    '''

    ML_obj = MLclassifier(bcftools_folder=pytest.config.getoption('bcftools_folder'))
    outfile = ML_obj.train(outprefix="data/outdir/fitted_logreg_indels",
                           tp_annotations=pytest.config.getoption('--tp_annotations_indels'),
                           fp_annotations=pytest.config.getoption('--fp_annotations_indels'))

    assert os.path.isfile(outfile) is True

def test_apply_model():

    ML_obj = MLclassifier(bcftools_folder=pytest.config.getoption('bcftools_folder'),
                          fitted_model='data/outdir/fitted_logreg_snps.sav')

    outfile = ML_obj.predict(outprefix="data/outdir/predictions",
                             annotation_f=pytest.config.getoption('--tp_annotations_snps'),
                             cutoff=0.95)

    assert os.path.isfile(outfile) is True

def test_rfe(clean_tmp):

    ML_obj = MLclassifier(bcftools_folder=pytest.config.getoption('bcftools_folder'))

    select_feats = ML_obj.rfe(tp_annotations=pytest.config.getoption('--tp_annotations_indels'),
                              fp_annotations=pytest.config.getoption('--fp_annotations_indels'),
                              n_features=5)

    assert all([a == b for a, b in zip(select_feats, ['[4]IDV', '[5]IMF', '[6]VDB',
                                                      '[7]SGB', '[11]HOB'])])
