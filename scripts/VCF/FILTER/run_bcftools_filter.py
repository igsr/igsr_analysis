from VcfFilter import VcfFilter

import argparse
import os

#get command line arguments

parser = argparse.ArgumentParser(description='Script to apply bcftools filter on a VCF file')

#parameters
parser.add_argument('--bcftools_folder', type=str, required=True, help='Folder containing the Bcftools binary' )
parser.add_argument('--filename', type=str, required=True, help='Name (without the fullpath) of the VCF file that will be analysed. It assumes that the filename format is for example lc_bams.gatk.xxxx.vcf.gz, where lc_bams is the analysis group and gatk is the method used' )
parser.add_argument('--filter_name', type=str, required=True, help='Filter name used for annotate the VCF' )
parser.add_argument('--expression', type=str, required=True, help='Expression used for the filtering' )

args = parser.parse_args()

if __name__ == '__main__':
    
    vcf_f=VcfFilter(vcf=args.filename,bcftools_folder=args.bcftools_folder)
    vcf_f.bcftools_filter(name=args.filter_name,expression=args.expression)
    
