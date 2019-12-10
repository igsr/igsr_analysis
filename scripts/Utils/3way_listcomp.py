import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--refids', required=True, help='File with ids that will be used as the reference and will be compared with', type=str)
parser.add_argument('--file_list', nargs='+', help='List with files containing the ids that will be compared with --refids', required=True)

args = parser.parse_args()

ref_ids=args.refids
filelist=args.file_list

def get_ids_from_file(file):
    ids=[]
    with open(file) as f:
        for line in f:
            line=line.rstrip()
            ids.append(line)
    return ids

refids=get_ids_from_file(ref_ids)

idlist=[]
for f in filelist:
    idlist.append(get_ids_from_file(f))
    
header="ref\t"+"\t".join(filelist)
print(header)

for id in refids:
    idstr=id+"\t"
    for idl in idlist:
        if id in idl:
            idstr+="1\t"
        else:
            idstr+="0\t"
    print(idstr)

