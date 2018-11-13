import gzip
import sys
import re
import argparse

parser = argparse.ArgumentParser(description='Script to replace the chromosome numbers by ucsc-type chr ids (i.e. chr1, chr2, ....')
parser.add_argument('--input', type=str, required=True, help='vcf.gz input file to recode' )
parser.add_argument('--output', type=str, required=True, help='vcf.gz that will be used to write the output' )
parser.add_argument('--ucsc', type=bool, required=False, help='Input chro names are in ucsc format and will be converted to ensembl-like format' )

args = parser.parse_args()

p = re.compile( '^chr' )

f=gzip.open(args.output,'wb');

with gzip.open(args.input,'r') as fin:
    for line in fin:
        if not line.startswith(b"#"):
            bits=line.split(b"\t")
            chrname=bits[0].decode("utf-8")
            nchrname=""
            if args.ucsc==True:
                m = p.match( chrname )
                if not m:
                    raise Exception("This file does not have UCSC-style chros")
                nchrname = re.sub('chr', '', chrname)
            else:                
                nchrname="chr"+chrname
            bits[0]=nchrname.encode('utf-8')
            nline=b'\t'.join(bits)
            f.write(nline)
        else:
            f.write(line)
f.close()
