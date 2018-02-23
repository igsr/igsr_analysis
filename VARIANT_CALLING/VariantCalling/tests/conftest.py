import pytest

def pytest_addoption(parser):
    parser.addoption('--bam', default='data/exampleBAM.bam', action='store_true', help='Path to test BAM file')
    parser.addoption('--reference', default='data/exampleFASTA.fasta', action='store_true', help='Path to Fasta file')
    parser.addoption('--gatk_folder', default='~/bin/GATK/', action='store_true', help='Path to folder containing the GATK jar file')
    parser.addoption('--bgzip_folder', default='/nfs/production/reseq-info/work/ernesto/bin/anaconda3/bin/', action='store_true', help='Path to folder containing the Bgzip binary')


