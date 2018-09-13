from VcfNormalize import VcfNormalize

import argparse
import os

#get command line arguments

parser = argparse.ArgumentParser(description='Script to run GATK VariantsToAllelicPrimitives in order to decompose MNPs into more basic/primitive alleles')

parser.add_argument('--gatk_folder', type=str, required=True, help='Folder containing GATK jar file' )
parser.add_argument('--bgzip_folder', type=str, required=True, help='Folder containing bgzip' )
parser.add_argument('--vcf', type=str, required=True, help='Path to the VCF file that will be analysed' )
parser.add_argument('--outprefix', type=str, required=True, help='Prefix for output file' )
parser.add_argument('--reference', type=str, required=True, help='Path to the reference Fasta file' )
parser.add_argument('--compress', type=str, required=False, help='Compress the output file' )

args = parser.parse_args()

if __name__ == '__main__':
        
    vcfallprim = VcfNormalize(vcf=args.vcf,gatk_folder=args.gatk_folder,bgzip_folder=args.bgzip_folder)
    vcfallprim.run_gatk_VariantsToAllelicPrimitives(outprefix=args.outprefix,reference=args.reference,compress=args.compress)
    
