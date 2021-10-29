import os
import pytest

from VCF.VCFfilter.MLclassifier import MLclassifier

def test_train_snps(datadir):
    """
    Train the model for SNPs
    """

    ML_obj = MLclassifier()
    outfile = ML_obj.train(outprefix="{0}/outdir/fitted_logreg_snps".format(datadir),
                           tp_annotations="{0}/TP_annotations_snps.chr20.tsv".format(datadir),
                           fp_annotations="{0}/FP_annotations_snps.chr20.tsv".format(datadir))

    assert os.path.isfile(outfile) is True

def test_train_snps_gz(datadir):
    """
    Train the model for SNPs using 2 gzipped annotation files
    """

    ML_obj = MLclassifier()
    outfile = ML_obj.train(outprefix="{0}/outdir/fitted_logreg_snps".format(datadir),
                           tp_annotations="{0}/TP_annotations_snps.chr20.tsv.gz".format(datadir),
                           fp_annotations="{0}/FP_annotations_snps.chr20.tsv.gz".format(datadir))

def test_train_indels(datadir):
    """
    Train the model for INDELs
    """

    ML_obj = MLclassifier()
    outfile = ML_obj.train(outprefix="{0}/outdir/fitted_logreg_indels".format(datadir),
                           tp_annotations="{0}/TP_annotations_indels.chr20.tsv".format(datadir),
                           fp_annotations="{0}/TP_annotations_indels.chr20.tsv".format(datadir))

    assert os.path.isfile(outfile) is True

def test_apply_model(clean_tmp, datadir):

    ML_obj = MLclassifier(fitted_model="{0}/outdir/fitted_logreg_snps.sav".format(datadir))

    outfile = ML_obj.predict(outprefix="{0}/outdir/predictions".format(datadir),
                             annotation_f="{0}/TP_annotations_snps.chr20.tsv".format(datadir),
                             cutoff=0.95)

    assert os.path.isfile(outfile) is True

def test_rfe(datadir, clean_tmp):

    ML_obj = MLclassifier()

    select_feats_report = ML_obj.rfe(tp_annotations="{0}/TP_annotations_indels.chr20.tsv".format(datadir),
        fp_annotations="{0}/TP_annotations_indels.chr20.tsv".format(datadir),
        n_features=5,
        outreport="{0}/outdir/out_rfe.txt".format(datadir))

    assert os.path.isfile(select_feats_report) is True