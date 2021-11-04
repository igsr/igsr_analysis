import pytest
import pdb
import os

from GRAPH.VgToolkit import VG

@pytest.fixture
def vg_object(datadir):
    '''Returns a VG object'''
    vg_object = VG()

    return vg_object

def test_run_autoindex(vg_object, datadir):
    '''
    Test function to run 'vg autoindex'
    '''

    ofiles = vg_object.run_autoindex(ref=f"{datadir}/VG/hs37d5.fa",
                                     vcf=f"{datadir}/VG/chr20.vcf.gz",
                                     prefix=f"{datadir}/outdir/test.autoindex")
    
    for f in [f"{os.getcwd()}/data/outdir/test.autoindex.dist", f"{os.getcwd()}/data/outdir/test.autoindex.giraffe.gbz", 
    f"{os.getcwd()}/data/outdir/test.autoindex.min"]:
        assert f in ofiles

def test_run_chunk(vg_object, datadir):
    '''
    Test function to run 'vg chunk'
    '''

    outfiles = vg_object.run_chunk('M',
                                  x=f"{datadir}/VG/x.xg",
                                  O="pg",
                                  b=f"{datadir}/outdir/test")

    e_ofiles = [f"{datadir}/outdir/test_x.pg"]
    assert all([a == b for a, b in zip(outfiles, e_ofiles)])

def test_run_snarls(vg_object, datadir):
    '''
    Test  function to run 'vg snarls'
    '''

    outfile = vg_object.run_snarls(graph_f = f"{datadir}/outdir/test_x.pg",
                                   prefix=f"{datadir}/outdir/test_x")
    
    assert os.path.isfile(outfile) is True

def read(vg_object, datadir):
    '''
    Test function to run 'vg giraffe' on a FASTQ file
    '''
    
    outfile = vg_object.run_giraffe(Z=f"{datadir}/outdir/test.autoindex.giraffe.gbz",
                                    m=f"{datadir}/outdir/test.autoindex.min", 
                                    d=f"{datadir}/outdir/test.autoindex.dist", 
                                    fastq=f"{datadir}/VG/input1.fq", 
                                    prefix=f"{datadir}/outdir/test")

    assert os.path.isfile(outfile) is True

def test_run_giraffe_paired(vg_object, datadir):
    '''
    Test function to run 'vg giraffe' on a pair of mated FASTQ files
    '''
    
    outfile = vg_object.run_giraffe(Z=f"{datadir}/outdir/test.autoindex.giraffe.gbz",
                                    m=f"{datadir}/outdir/test.autoindex.min", 
                                    d=f"{datadir}/outdir/test.autoindex.dist", 
                                    fastq=f"{datadir}/VG/input1.fq, {datadir}/VG/input2.fq",
                                    prefix=f"{datadir}/outdir/test")

    assert os.path.isfile(outfile) is True

def test_run_stats(vg_object, datadir):
    '''
    Test function to run 'vg stats' on a .gam file
    '''
    outfile = vg_object.run_stats(aln_f=f"{datadir}/outdir/test.gam")

    assert os.path.isfile(outfile) is True

def test_run_augment(vg_object, datadir):
    '''
    Test function to run 'vg augment' on a .gam file
    '''
    outfiles = vg_object.run_augment(graph_f=f"{datadir}/VG/x.vg",
                                     aln_f=f"{datadir}/VG/aln.gam",
                                     prefix=f"{datadir}/outdir/aln")

    assert os.path.isfile(outfiles[0]) is True
    assert os.path.isfile(outfiles[1]) is True

def test_run_pack(vg_object, datadir):
    '''
    Test function to run 'vg pack' on a .gam file
    '''
    outfile = vg_object.run_pack(vg_f=f"{datadir}/VG/x.vg",
                                 aln_f=f"{datadir}/VG/aln.gam",
                                 prefix=f"{datadir}/outdir/aln",
                                 Q=5)
    
    assert os.path.isfile(outfile) is True

def test_run_call(vg_object, datadir, clean_tmp):
    '''
    Test function to run 'vg call' to generate a VCF file
    '''
    outfile = vg_object.run_call(graph_f=f"{datadir}/VG/x.vg",
                                 k=f"{datadir}/outdir/aln.pack",
                                 prefix=f"{datadir}/outdir/aln")
    
    assert os.path.isfile(outfile) is True