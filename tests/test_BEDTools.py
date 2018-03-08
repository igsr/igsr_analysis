import os
import pytest
import glob
import warnings
from BEDTools import BEDTools

# test_BEDTools.py

@pytest.fixture
def bedtools_object(scope='module'):
    '''Returns a BEDTools  object'''

    print("Creating the object\n")
    bedtools_folder = pytest.config.getoption("--bedtools_folder")
    bedtools_object=BEDTools(bedtools_folder=bedtools_folder)
    return bedtools_object

def test_make_windows(bedtools_object):
    coordlist=bedtools_object.make_windows(g='data/chr1.genome', w=100000000)
    
    assert coordlist[0]==['chr1', '0', '100000000']
    assert coordlist[1]==['chr1', '100000000', '200000000']
    assert coordlist[2]==['chr1', '200000000', '248956422']

def test_make_windows_and_subtract(bedtools_object):
    coordlist=bedtools_object.make_windows(g='data/chr1.genome', subtract='data/subtract.bed',w=100000000)

    assert coordlist[0]==['chr1', '0', '100000000']
    assert coordlist[1]==['chr1', '100000000', '100000100']
    assert coordlist[2]==['chr1', '100000200', '100000500']
    assert coordlist[3]==['chr1', '100000600', '200000000']
    assert coordlist[4]==['chr1', '200000000', '248956422']

