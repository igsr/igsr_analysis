# This script will take a list of files with the following format:
# chr1\t1000 and will generate a matrix of presence/absence of each of the
# sites in each of the files. 1 will mean presence and 0 will mean absence
# Example of the output:
# #chr pos file1.txt file2.txt
# chr1 1 1 1
# chr1 10 1 0
# chr1 11 0 1
# chr1 20 1 1
#
# Usage: python generate_matrix_pos.py --file_list file1.txt file2.txt

import argparse
import pdb
import os

parser = argparse.ArgumentParser()

parser.add_argument('--file_list', nargs='+', help='List with files containing the ids that will be compared with --refids', required=True)

args = parser.parse_args()

filelist=args.file_list

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

set_dict=AutoVivification()
all_pos=AutoVivification()

def parse_file(path):

    #get basename that will be used as id for this dataset
    dataset_id='.'.join(os.path.basename(path).split('.')[:-1])
    with open(path) as f:
        for line in f:
            line=line.rstrip()
            (chr,pos)=line.split("\t")
            pos=int(pos)
            set_dict[dataset_id][chr][pos]=0
            all_pos[chr][pos]=0
            
            
for f in filelist:
    parse_file(path=f)

datasets=set_dict.keys()

header="#chr pos "
header+=" ".join(datasets)
print(header)

chr_order=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']

for c in chr_order:
    sorted_pos=sorted(all_pos[c])
    for pos in sorted_pos:
        line="{0} {1} ".format(c,pos)
        for dataset_id in datasets:
            if pos in set_dict[dataset_id][c]:
                line+="1 "
            else:
                line+="0 "
        print(line)
 
