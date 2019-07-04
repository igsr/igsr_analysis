import glob
import os
import pdb
from collections import defaultdict
from prettytable import PrettyTable

all_blocks_d = defaultdict(list)
larg_block_d = defaultdict(list)

for file in glob.glob("/nfs/production/reseq-info/work/ernesto/isgr/VARIANT_CALLING/VARCALL_ALLGENOME_13022017/BENCHMARKING/WHATSHAP/SNPS/GRCh38/out.*.log"):
    chromosome=int(os.path.basename(file).split('.')[1])
    all_block=False
    larg_block=False
    with open(file) as f:
        for line in f:
            line=line.rstrip("\n")
            line=line.lstrip()
            elms=line.split(":")
            if len(elms)>1:
                if elms[0]=="ALL INTERSECTION BLOCKS":
                    all_block=True
                elif all_block is True and elms[0]!="LARGEST INTERSECTION BLOCK":
                    key=elms[0].strip()
                    value=elms[1].strip()
                    all_blocks_d[chromosome].append((key,value))
                elif elms[0]=="LARGEST INTERSECTION BLOCK":
                    larg_block=True
                    all_block=False
                elif larg_block is True and elms[0]!="Hamming distance [%]":
                    key=elms[0].strip()
                    value=elms[1].strip()
                    larg_block_d[chromosome].append((key,value))
                elif larg_block is True and elms[0]=="Hamming distance [%]":
                    key=elms[0].strip()
                    value=elms[1].strip()
                    larg_block_d[chromosome].append((key,value))
                    break

all_blocks_tbl = PrettyTable()
all_blocks_tbl.field_names = ["chr", "Phased pairs of variants assessed", "Switch errors", 
                              "Switch error rate", 'Switch/flip decomposition', 
                              'Switch/flip rate', 'Block-wise Hamming distance', 
                              'Block-wise Hamming distance [%]']

for chrom in sorted(all_blocks_d.keys()):
    values= all_blocks_d[chrom]
    vs=[chrom,values[0][1],values[1][1],values[2][1],values[3][1],
        values[4][1],values[5][1],values[6][1]]
    all_blocks_tbl.add_row(vs)

print("#\n## All intersection blocks##\n#")
print(all_blocks_tbl)    

larg_block_tbl = PrettyTable()
larg_block_tbl.field_names = ["chr", "Phased pairs of variants assessed", "Switch errors",
                              "Switch error rate", 'Switch/flip decomposition',
                              'Switch/flip rate', 'Hamming distance',
                              'Hamming distance [%]']

for chrom in sorted(all_blocks_d.keys()):
    values= all_blocks_d[chrom]
    vs=[chrom,values[0][1],values[1][1],values[2][1],values[3][1],
        values[4][1],values[5][1],values[6][1]]
    larg_block_tbl.add_row(vs)

print("#\n## Largest intersection block##\n#")
print(larg_block_tbl)
