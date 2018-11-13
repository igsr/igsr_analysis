from VcfQC import VcfQC

import argparse
import os
import logging

#get command line arguments

parser = argparse.ArgumentParser(description='Script to run bcftools stats on a certain VCF file')

parser.add_argument('--bcftools_folder', type=str, required=True, help='Folder containing the VCFtools binary' )
parser.add_argument('--filename', type=str, required=True, help='Path to the VCF file that will be analysed' )
parser.add_argument('--outpath', type=str, required=True, help='Path for output file')

args = parser.parse_args()

if __name__ == '__main__':
    
    vcfQC = VcfQC(vcf=args.filename,bcftools_folder=args.bcftools_folder)
    stats=vcfQC.stats(filter_str="PASS,.",outpath=args.outpath,verbose=True)


