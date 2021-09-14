import pytest

@pytest.fixture
def vg_object(datadir, bcftools_folder):
    '''Returns a VG object'''

    vg_object = VG()

    return vg_object

def test_run_giraffe(vg_object):
    '''
    Test function to run 'vg giraffe' on a FASTQ file
    '''
    pdb.set_trace()

    vg_object.run_giraffe()