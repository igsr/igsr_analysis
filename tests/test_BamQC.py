import os
import pytest
import glob
import warnings
from BamQC import BamQC
import pdb

# test_BamQC.py
@pytest.fixture
def bam_object(scope='module'):
    '''Returns an  object'''
    print("Creating the object\n")
    bam_file = pytest.config.getoption("--bam")
    chk_indel_folder = pytest.config.getoption("--chk_indel_folder")
    samtools_folder = pytest.config.getoption("--samtools_folder")
    java_folder = pytest.config.getoption("--java_folder")
    picard_folder = pytest.config.getoption("--picard_folder")

    bam_object = BamQC(bam=bam_file,chk_indel_folder=chk_indel_folder,
                       samtools_folder=samtools_folder,java_folder=java_folder,
                       picard_folder=picard_folder)

    return bam_object

@pytest.fixture
def clean_tmp():
    yield
    print("Cleanup files")
    files = glob.glob('data/outdir/*')
    for f in files:
        os.remove(f)

def test_get_simple_stats(bam_object):
    stats=bam_object.get_simple_stats()
    assert stats['total_no_reads']==33
    assert stats['no_duplicates']==0
    assert stats['total_no_mapped']==33
    assert stats['no_properly_paired']==32

def test_get_contigs(bam_object):
    contigs=bam_object.get_contigs()
    assert len(contigs.keys())==1

def test_list_of_samples(bam_object):
    samples=bam_object.list_of_samples()
    assert samples[0]=='exampleBAM.bam'

def test_list_of_readgroups(bam_object):
    rgroups=bam_object.list_of_readgroups()
    assert rgroups[0]=='exampleBAM.bam'

def test_run_samtools_depth(bam_object):
    depthO_l=bam_object.run_samtools_depth(chros='chr1')
    assert depthO_l[0].contig=='chr1'
    assert depthO_l[0].bases_mapped==2052
    assert depthO_l[0].breadth==0.02052
    assert depthO_l[0].depth==0.02508
    assert depthO_l[0].sum_of_depths==2508
    assert depthO_l[0].length==100000
    assert depthO_l[0].max==3

def test_aggregate_stats(bam_object):
    depthO_l=bam_object.run_samtools_depth(chros='chr1')
    depthO=bam_object.aggregate_stats(depthO_l)
    assert depthO.contig=='chr1'
    assert depthO.bases_mapped==2052
    assert depthO.breadth==0.02052
    assert depthO.depth==0.02508
    assert depthO.sum_of_depths==2508
    assert depthO.length==100000
    assert depthO.max==3

def test_run_CollectHsMetrics(bam_object):
    cMetrics=bam_object.run_CollectHsMetrics(baits_file='data/test.ival')
    assert cMetrics.metrics['TOTAL_READS']=='33'
    assert cMetrics.metrics['PF_UNIQUE_READS']=='33'
    assert cMetrics.metrics['PF_UQ_BASES_ALIGNED']=='2508'

def test_Wes_create_cov_barplot(bam_object):
    cMetrics=bam_object.run_CollectHsMetrics(baits_file='data/test.ival',outfile='data/outdir/example')
    cMetrics.create_cov_barplot('data/outdir/example.CollectHsMetrics.barplot.pdf')

def test_run_CollectWgsMetrics(bam_object):
    cMetrics=bam_object.run_CollectWgsMetrics(reference='data/exampleFASTA.fasta')
    assert cMetrics.metrics['MEAN_COVERAGE']=='0.00139'
    assert cMetrics.metrics['PCT_1X']=='0.00139'
    assert cMetrics.metrics['SD_COVERAGE']=='0.037257'

def test_Wgs_create_cov_barplot(bam_object, clean_tmp):
    cMetrics=bam_object.run_CollectWgsMetrics(reference='data/exampleFASTA.fasta',outfile='data/outdir/out_runCollectWgsMetrics.txt')
    cMetrics.create_cov_barplot('data/outdir/example.CollectWgsMetrics.barplot.pdf')


