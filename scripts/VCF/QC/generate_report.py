import argparse
import glob
import re
import pdb
import os
import pandas as pd
import openpyxl

from tabulate import tabulate

def check_variantype(value):
    if value !='snps' and value!='indels':
        raise argparse.ArgumentTypeError("%s is an invalid variant type" % value)
    return value


parser = argparse.ArgumentParser(description='Script to generate report on the benchmarking of a VCF')

parser.add_argument('--dir', required=True, help='Folder containing the per-chr folders generated by compare_with_giab.nf')
parser.add_argument('--vt', required=True, help='Type of variant to analyze. Possible values are \'snps\' and \'indels\'', type=check_variantype)
parser.add_argument('--outfile', required=True, help='Name for output spreadsheet. Example: out.xlsx')
parser.add_argument('--subset', required=True, help='Generate report with High conf sites or with all sites. Possible values are: \'highconf\' or \'all\'')


args = parser.parse_args()


p=re.compile(".*/results_(.*)")


class BcftoolsStats(object):
    '''
    Class to store the results of running BCFtools stats on a VCF file
    '''

    def __init__(self, filename=None, summary_numbers=None, ts_tv=None,
                 ts_tv_1stalt=None, no_singleton_snps=None):
        '''
        Constructor
        Parameters
        ----------
        filename : str
                   Filename of the VCF that was used to run bcftools stats
        summary_numbers : dict
                          Dictionary containing the basic stats. i.e.:
                              number of samples:      1
                              number of records:      1867316
                              .....
        ts_tv : float
                ts/tv ratio
        ts_tv_1stalt : float
                       ts/tv (1st ALT)
        no_singleton_snps : int
        '''

        self.filename = filename
        self.summary_numbers = summary_numbers
        self.ts_tv = ts_tv
        self.ts_tv_1stalt = ts_tv_1stalt
        self.no_singleton_snps = no_singleton_snps

    def __str__(self):
        sb = []
        for key in self.__dict__:
            sb.append("{key}='{value}'".format(key=key, value=self.__dict__[key]))

        return ', '.join(sb)

    def __repr__(self):
        return self.__str__()



def parse_stats_file(f):
    '''
    Function to parse a stats file
    :param f:

    Returns
    -------
    BcftoolsStats object
    '''

    stats = BcftoolsStats(filename=f)

    with open(f) as fi:
        d = {}
        for line in fi:
            line = line.rstrip('\n')
            if line.startswith('SN\t'):
                key = line.split('\t')[2]
                value = int(line.split('\t')[3])
                d[key] = value
            elif line.startswith('TSTV\t'):
                ts_tv = line.split('\t')[4]
                ts_tv_1stalt = line.split('\t')[7]
                stats.ts_tv = ts_tv
                stats.ts_tv_1stalt = ts_tv_1stalt
            elif line.startswith('SiS\t'):
                no_singleton_snps = line.split('\t')[3]
                stats.no_singleton_snps = no_singleton_snps

    stats.summary_numbers = d
    return stats

data = dict()

for dir in glob.glob(args.dir+"/results_*"):
    print("Processing: {0}".format(dir))
    m = p.match(dir)
    if m:
        chr=m.group(1)
        numbers=dict()
        res=None
        if args.subset=="highconf":
            res = [f for f in glob.glob(dir+"/*.stats") if "highconf" in f]
        elif args.subset=="all":
            res = [f for f in glob.glob(dir+"/*.stats") if not "highconf" in f]
        else:
            raise Exception("Value not recognised for subset argument: {0}".format(args.subset))
        for f in res:
            type=os.path.basename(f).split(".")[0]
            bcfobj=parse_stats_file(f)
            sum_dict=bcfobj.summary_numbers
            if args.vt=='snps':
                numbers[type]=sum_dict['number of SNPs:']
            elif args.vt=='indels':
                numbers[type] = sum_dict['number of indels:']
        chr_stripped=int(chr.replace("chr",""))
        if chr_stripped is not 'X': int(chr_stripped)
        data[chr_stripped]=numbers

    else:
        raise Exception('No chromosome was fetched from dir name')


df = pd.DataFrame.from_dict(data, orient='index')
df['total_cat1']=df.TP+df.FN
df['total_cat2']=df.TP+df.FP
df['%_TP']=round(df.TP*100/df.total_cat1,2)
df['%_FN']=round(df.FN*100/df.total_cat1,2)
df['%_FP']=round(df.FP*100/df.total_cat2,2)

df1=df.sort_index()

print(tabulate(df1, headers='keys', tablefmt='psql'))


writer = pd.ExcelWriter(args.outfile)
df1.to_excel(writer, 'Benchmarking')
writer.save()
