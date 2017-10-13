import pytest

def pytest_addoption(parser):
    parser.addoption('--hostname', default='mysql-g1kdcc-public', action='store_true', help='host name')
    parser.addoption('--username', default='g1krw', action='store_true', help='user name')
    parser.addoption('--port', default=4197, action='store_true', help='port')
    parser.addoption('--pwd', default='test', action='store_true', help='password')
    parser.addoption('--db', default='g1k_archive_staging_track', action='store_true', help='database name')

