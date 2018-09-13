from VcfQC import VcfQC

import argparse
import os
import logging

#get command line arguments

parser = argparse.ArgumentParser(description='Script to get the chros that are present in a VCF file. If filters are present in the VCF, then only the chros for the filtered variants will be considered')

parser.add_argument('--bcftools_folder', type=str, required=True, help='Folder containing the VCFtools binary' )
parser.add_argument('--filename', type=str, required=True, help='Path to the VCF file that will be analysed' )
parser.add_argument('--filter_str', type=str, required=False, help='i.e. "PASS,." If defined, then apply the filters before getting the list of chros' )

args = parser.parse_args()

if __name__ == '__main__':
    
    vcfQC = VcfQC(vcf=args.filename,bcftools_folder=args.bcftools_folder)

    if args.filter_str:
        chros_list=vcfQC.get_chros(filter_str=args.filter_str,verbose=True)
    else:
        chros_list=vcfQC.get_chros(verbose=True)

