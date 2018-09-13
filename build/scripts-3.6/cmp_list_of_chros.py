from VcfQC import VcfQC

import argparse
import os
import logging

#get command line arguments

parser = argparse.ArgumentParser(description='Script to get the chros that are present in a VCF file (after applying the filters specified using the --filter_str arg) and compare the obtained list with a list of chros passed in a file')

parser.add_argument('--bcftools_folder', type=str, required=True, help='Folder containing the VCFtools binary' )
parser.add_argument('--filename', type=str, required=True, help='Path to the VCF file that will be analysed' )
parser.add_argument('--chr_f', type=str, required=True, help='File with chromosomes (one chromosome per line) for which the list of chros in the VCF will be compared' )
parser.add_argument('--filter_str', type=str, required=False, help='i.e. "PASS,." If defined, then apply the filters before getting the list of chros' )

args = parser.parse_args()

if __name__ == '__main__':
    
    vcfQC = VcfQC(vcf=args.filename,bcftools_folder=args.bcftools_folder)

    chros_dict={}
    if args.filter_str:
        chros_dict=vcfQC.get_chros(filter_str=args.filter_str,chr_f=args.chr_f,verbose=True)
    else:
        chros_dict=vcfQC.get_chros(verbose=True,chr_f=args.chr_f)
    if len(chros_dict['in_A'])>0:
        print("FAILED!\n")
        print("Chros present in the VCF and not in chr_file:")
        print(','.join(chros_dict['in_A']))
    
    if len(chros_dict['in_B'])>0:
        print("FAILED!\n")
        print("Chros not in the VCF and present in chr_file:")
        print(','.join(chros_dict['in_B']))
