from VcfQC import VcfQC

import argparse
import os
import logging

#get command line arguments

parser = argparse.ArgumentParser(description='Script to get the number of variants falling within a particular region')

parser.add_argument('--bedtools_folder', type=str, required=True, help='Folder containing the BEDtools binary' )
parser.add_argument('--filename', type=str, required=True, help='Path to the VCF file that will be analysed' )
parser.add_argument('--region', type=str, required=True, help='BED file with regions to be analyzed' )
parser.add_argument('--outprefix', type=str, required=True, help='Prefix for output file' )

args = parser.parse_args()

if __name__ == '__main__':
    
    vcfQC = VcfQC(vcf=args.filename,bedtools_folder=args.bedtools_folder)
    vcfQC.number_variants_in_region(region=args.region,outprefix)
   
