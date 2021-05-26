import os
import pytest
import glob
import shutil

from VCF.VCFIntegration import SNPTools

# test_vcfqc.py
@pytest.fixture
def vcf_object(datadir):
    """Returns a SNPTools object"""
    
    vcf_file = "{0}/test.vcf.gz".format(datadir)

    snptools_folder = None # folder containing the snptools binary
    if shutil.which('bamodel') is None:
        raise Exception("'bamodel' needs to by in $PATH")
    else:
        snptools_folder = os.path.dirname(shutil.which('bamodel'))

    vcf_object = SNPTools(vcf=vcf_file, snptools_folder=snptools_folder)
    
    return vcf_object

def test_run_bamodel(vcf_object, datadir, clean_tmp):
    
    vcf_object.run_bamodel(sample="NA12878_chr1_1000000_1001000",
                           bamfiles="{0}/SNPTools/bamlist.txt".format(datadir),
                           outdir="{0}/SNPTools/outdir/".format(datadir))

    assert os.path.exists("{0}/SNPTools/outdir/NA12878_chr1_1000000_1001000.raw".format(datadir))

def test_run_poprob(vcf_object, datadir, clean_tmp):
    
    vcf_object.run_poprob(outprefix="NA12878_chr1_1000000_1001000",
                          rawlist="{0}/SNPTools/rawlist.txt".format(datadir),
                          outdir="{0}/SNPTools/outdir/".format(datadir))
    assert os.path.exists("{0}/SNPTools/outdir/NA12878_chr1_1000000_1001000.prob".format(datadir))

def test_run_prob2vcf(vcf_object, datadir, clean_tmp):
    vcf_object.run_prob2vcf(probf="{0}/SNPTools/outdir/NA12878_chr1_1000000_1001000.prob".format(datadir),
                            outprefix="NA12878_chr1_1000000_1001000",
                            outdir="{0}/SNPTools/outdir/".format(datadir),
                            chro="chr1")
    assert os.path.exists("{0}/SNPTools/outdir/NA12878_chr1_1000000_1001000.vcf.gz".format(datadir))
