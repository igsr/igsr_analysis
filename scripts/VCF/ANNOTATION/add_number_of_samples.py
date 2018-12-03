import subprocess
import argparse

parser = argparse.ArgumentParser(description='Script to add the NS (number of samples) annotation to the annotation table')

parser.add_argument('--outfile',help='Output file name',required=True)
parser.add_argument('--file1',help='Path to table that will be used in order to add the annotation as the last column',required=True)
parser.add_argument('--phased_vcf',help='Path to phased VCF that will be used to obtain the number of samples',required=True)

args = parser.parse_args()

def add_number_of_samples(file1,ann_vcf):
    '''
    Function to add the NS (number of samples) annotation to the annotation table

    Args
    ----
    file1: string
           Path to table that will be used in order to add the annotation as the last column
    phased_vcf: string
                Path to phased VCF that will be used to obtain the number of samples

    Returns
    -------
    Path to the file that will contain the additional column

    '''

    cmd="bcftools query -l {0} |wc -l".format(args.phased_vcf)

    try:
        number_of_samples = subprocess.check_output(cmd, shell=True).decode("utf-8").rstrip('\n')
        wf=open(args.outfile,'w')
        with open(file1) as f:
            for line in f:
                line=line.rstrip("\n")
                if line.startswith('#'):
                    line=line+"\tNS"
                    wf.write(line+"\n")
                    continue
                else:
                    line= line+"\t{0}".format(number_of_samples)
                    wf.write(line+"\n")
    except subprocess.CalledProcessError as e:
        raise

add_number_of_samples(args.file1, args.phased_vcf)

