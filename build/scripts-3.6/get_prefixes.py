import argparse
import re

#get command line arguments

parser = argparse.ArgumentParser(description='Script to modify a ped file so the pedigree info is added')

#parameters
parser.add_argument('--ped', type=str, required=True, help='PED file that will be modified' )
parser.add_argument('--ped_info', type=str, required=True, help='PED file containing the pedigree info' )
parser.add_argument('--outf', type=str, required=True, help='Output file' )
args = parser.parse_args()

outf=open(args.outf,'w');

noped_counter=0
with open(args.ped) as f:
    for line in f:
        line=line.rstrip('\n')
        e=line.split("\t")
        sample=e[1]
        prefix=""
        with open(args.ped_info) as f1:
            seen=0
            for line1 in f1:
                e1=line1.split("\t")
                if e1[1]==sample:
                    seen=1
                    prefix="{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t".format(e1[0],e1[1],e1[2],e1[3],e1[4],e1[5])
                    
            if seen==0:
                noped_counter+=1
                prefix="{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t".format(e[0],e[1],e[2],e[3],e[4],e[5])
                        
        newline=prefix+"\t".join(e[6:])
        outf.write(newline+"\n")

outf.close

print("Number of variants without ped info: {0}".format(noped_counter))
