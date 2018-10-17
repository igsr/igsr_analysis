import argparse

parser = argparse.ArgumentParser(description='Script to add the Variant Type annotation to the VCF')

parser.add_argument('--outfile',help='Output file',required=True)
parser.add_argument('--file1',help='Path to table that will be used in order to add the annotation as the last column',required=True)
parser.add_argument('--header',help='Name used in the header for this annotation',required=True)
parser.add_argument('--label',help='Label that will be used and will be added as the last column in the table',required=True)

args = parser.parse_args()

def add_annotation(file1,header,label):
    '''
    Function to add a particular label to all rows in the table

    Args
    ----
    file1: string
           Path to table that will be used in order to add the annotation as the last column
    header: string
            Name used in the header for this annotation
    label: string
           Label that will be used and will be added as the last column in the table

    Returns
    -------
    Nothing
    '''

    wf=open(args.outfile,'w')
    with open(file1) as f:
        for line in f:
            line=line.rstrip("\n")
            if line.startswith('#'):
                line=line+"\t{0}".format(header)
                wf.write(line+"\n")
                continue
            line=line.rstrip("\n")
            line= line+"\t{0}".format(label)
            wf.write(line+"\n")
    wf.close()

add_annotation(args.file1,args.header,args.label)
