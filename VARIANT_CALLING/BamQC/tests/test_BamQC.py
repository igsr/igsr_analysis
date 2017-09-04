import os
import pytest
import glob
import warnings
from BamQC import BamQC

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
    files = glob.glob('data/BamQC/*')
    for f in files:
        os.remove(f)
'''
TODO: Binary does not work in MAC, check in Linux

def test_chkindel_rg(bam_object):
    bam_object.run_chk_indel_rg(outfile="data/BamQC/test.chkindel_rg.txt")
    assert 0
'''

def test_get_simple_stats(bam_object):
    stats=bam_object.get_simple_stats()
    assert stats['total_no_reads']==672
    assert stats['no_duplicates']==14
    assert stats['total_no_mapped']==667
    assert stats['no_properly_paired']==645

def test_get_contigs(bam_object):
    contigs=bam_object.get_contigs()
    assert len(contigs.keys())==1

def test_list_of_samples(bam_object):
    samples=bam_object.list_of_samples()
    assert samples[0]=='NA19006'

def test_list_of_readgroups(bam_object):
    rgroups=bam_object.list_of_readgroups()
    assert rgroups[0]=='ERR251631' or rgroups[0]=='ERR251632' 
    assert rgroups[1]=='ERR251631' or rgroups[1]=='ERR251632'

def test_run_samtools_depth(bam_object):
    depthO_l=bam_object.run_samtools_depth(chros='chr1')
    assert depthO_l[0].contig=='chr1'
    assert depthO_l[0].bases_mapped==9139
    assert depthO_l[0].breadth==3.6709235803525486e-05
    assert depthO_l[0].depth==0.0002573984614865649
    assert depthO_l[0].sum_of_depths==64081
    assert depthO_l[0].length==248956422
    assert depthO_l[0].max==16

def test_aggregate_stats(bam_object):
    depthO_l=bam_object.run_samtools_depth(chros='chr1')
    depthO=bam_object.aggregate_stats(depthO_l)
    assert depthO.contig=='chr1'
    assert depthO.bases_mapped==9139
    assert depthO.breadth==3.6709235803525486e-05
    assert depthO.depth==0.0002573984614865649
    assert depthO.sum_of_depths==64081
    assert depthO.length==248956422
    assert depthO.max==16


def test_run_CollectHsMetrics(bam_object):
    cMetrics=bam_object.run_CollectHsMetrics(baits_file='data/test.ival')
    assert cMetrics.metrics['TOTAL_READS']=='672'
    assert cMetrics.metrics['PF_UNIQUE_READS']=='658'
    assert cMetrics.metrics['PF_UQ_BASES_ALIGNED']=='64081'

def test_run_run_CollectWgsMetrics(bam_object):
    cMetrics=bam_object.run_CollectWgsMetrics(reference='data/GRCh38_full_analysis_set_plus_decoy_hla.chr1.fa')
    assert 0
