'''
Script for calculating descriptive statistics on a file containing posisions and the aggrgated coverage for a set of samples.
This script uses the Dask library in order to cope with a huge amount of data

USAGE: python calc_median.py /path/to/out.cov
'''

import dask.dataframe as dd
import sys

df = dd.read_csv(sys.argv[1] , names=['chr','pos','cov'], sep='\t', blocksize=34000000)

print("Descriptors: {0}".format(df['cov'].describe().compute()))

