import dask.dataframe as dd

df = dd.read_csv('/hps/nobackup/production/reseq-info/ernesto/HIGHCOV/CALC_COVERAGE/PIPELINE/test/results/files/out.cov', names=['chr','pos','cov'], sep='\t', blocksize=34000000)

print("Descriptors: {0}".format(df['cov'].describe().compute()))

