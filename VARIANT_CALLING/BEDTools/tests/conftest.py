import pytest

def pytest_addoption(parser):
    parser.addoption('--bedtools_folder', default='/homes/ernesto/bin/bedtools-2.25.0/bin/', action='store_true', help='Folder with bedtools binary')
