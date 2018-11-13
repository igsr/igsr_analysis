from VcfNormalize import VcfNormalize

import argparse
import os

#get command line arguments

parser = argparse.ArgumentParser(description='Script to run vt normalize on a VCF file')

parser.add_argument('--vcflib_folder', type=str, required=True, help='Path to folder containing the different vcflib binaries' )
parser.add_argument('--bgzip_folder', type=str, required=False, help='Folder containing the Bgzip binary' )
parser.add_argument('--filename', type=str, required=True, help='Path to the file to be analyzed' )
parser.add_argument('--outprefix', type=str, required=True, help='Prefix for output file' )
parser.add_argument('--compress', type=bool, required=False, help='BGZIP compress the VCF output' )


args = parser.parse_args()

if __name__ == '__main__':
            
    vcfallprim = VcfNormalize(vcf=args.filename,vcflib_folder=args.vcflib_folder, bgzip_folder=args.bgzip_folder)
    vcfallprim.run_vcfallelicprimitives(outprefix=args.outprefix,compress=args.compress)
    
