import os
import pytest
from ReseqTrackDB import *

# test_ReseqTrackDB.py

@pytest.fixture
def reseqtrackdb_object():
    '''Returns an object'''
    hostname = pytest.config.getoption("--hostname")
    username = pytest.config.getoption("--username")
    port = pytest.config.getoption("--port")
    pwd = pytest.config.getoption("--pwd")
    db = pytest.config.getoption("--db")
    reseqtrackdb_object=ReseqTrackDB(host=hostname,user=username,port=port,pwd=pwd,db=db)
    return reseqtrackdb_object

def test_store(reseqtrackdb_object):
    
    #Create a File object with the new file to be stored
    f_obj = File(path='data/test1.vcf.gz', type='TEST_TYPE')
    res=f_obj.store(reseqtrackdb_object, dry=True, do_md5=True)

    assert res=='data/test1.vcf.gz'

def test_rename():
    '''
    This test will rename the File object with name= 'lc.bcftools.vcf.gz'
    to name='lc.bcftools.20171010.sites.vcf.gz'
    '''

    f_obj = File(path='data/lc.bcftools.vcf.gz', type='TEST_TYPE')

    f_obj.rename(filelayout=['set','caller','extension','compression'], newlayout=['set','caller'],
                 extension='sites.vcf', add_date=False, compression="gz")

    assert f_obj.name=='lc.bcftools.sites.vcf.gz'
