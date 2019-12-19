from Bio import SeqIO
import argparse
import pdb

parser = argparse.ArgumentParser(description='Script to generate a .BED file for a particular character type present in --fasta. This script can be used for example to detect \
the genomic coordinates for the stretches of characters produced after running the genome accessibility mask tool (bamUtil). It can also be used to get the coordinates of the \
ambiguous IUPAC bases that are present in --fasta')

parser.add_argument('--fasta', required=True, help='Fasta file to analyze')
parser.add_argument('--type', required=True, help='Character to search for. Examples are: N,L,H,Z,Q,R,M')

f=open('selected.bed','w');

args = parser.parse_args()

type=args.type

for seq_record in SeqIO.parse(args.fasta, "fasta"):
    coords=0 # starting by 0 to comply with .BED
    nts=str(seq_record.seq)
    nt_seen=False
    start=None
    end=None
    for nt in nts:
        if nt==type and nt_seen is False:
            start=coords
            nt_seen=True
        elif nt!=type and nt_seen is True:
            end=coords
            f.write("{0}\t{1}\t{2}\n".format(seq_record.id,start,end))
            nt_seen=False
        coords+=1
    if nt_seen is True:
        end=coords
        f.write("{0}\t{1}\t{2}\n".format(seq_record.id,start,end))

f.close
