from VcfQC import VcfQC

import argparse
import os

#get command line arguments

parser = argparse.ArgumentParser(description='Script to get the variant density for a certain VCF file')

parser.add_argument('--bedtools_folder', type=str, required=True, help='Folder containing the bedtools binary' )
parser.add_argument('--outprefix', type=str, required=True, help='Prefix for output' )
parser.add_argument('--filename', type=str, required=True, help='Path to the VCF file that will be analysed' )
parser.add_argument('--length', type=str, required=True, help='Window length (in bp) that will be used for calculating the density' )
parser.add_argument('--genome_file', type=str, required=True, help='File with genome sizes' )
parser.add_argument('--r_folder', type=str, required=False, help='Path to R binary' )
parser.add_argument('--r_scripts', type=str, required=True, help='Path to folder containing the R scripts required for constructing some of the plots used by this class' )

args = parser.parse_args()

if __name__ == '__main__':
        
    vcfQC = VcfQC(vcf=args.filename,bedtools_folder=args.bedtools_folder,r_folder=args.r_folder,r_scripts=args.r_scripts)
    
    vcfQC.plot_variant_density(length=args.length,genome=args.genome_file,outprefix=args.outprefix)
    
    

