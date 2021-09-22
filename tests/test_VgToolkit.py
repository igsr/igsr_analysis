import pytest
import pdb
import os

from GRAPH.VgToolkit import VG

@pytest.fixture
def vg_object(datadir, bcftools_folder):
    '''Returns a VG object'''
    VG.vg_folder="/nfs/production/flicek/ensembl/variation/software//vg-v1.34.0/"
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

def test_run_giraffe(vg_object, datadir):
    '''
    Test function to run 'vg giraffe' on a FASTQ file
    '''
    
    outfile = vg_object.run_giraffe(gbz_f=f"{datadir}/outdir/test.autoindex.giraffe.gbz",
                                    min=f"{datadir}/outdir/test.autoindex.min", 
                                    dist=f"{datadir}/outdir/test.autoindex.dist", 
                                    fastq=f"{datadir}/VG/input.fq", 
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
    outfiles = vg_object.run_augment(vg_f=f"{datadir}/VG/x.vg",
                                    aln_f=f"{datadir}/VG/aln.gam",
                                    prefix=f"{datadir}/outdir/aug")

    assert os.path.isfile(outfiles[0]) is True
    assert os.path.isfile(outfiles[1]) is True

