from VcfNormalize import VcfNormalize

import argparse
import os

#get command line arguments

parser = argparse.ArgumentParser(description='Script to run vt normalize on a VCF file')

parser.add_argument('--vt_folder', type=str, required=True, help='Folder containing the vt binary' )
parser.add_argument('--bgzip_folder', type=str, required=False, help='Folder containing the Bgzip binary' )
parser.add_argument('--filename', type=str, required=True, help='Path to the file to be analyzed' )
parser.add_argument('--outprefix', type=str, required=True, help='Prefix for output file' )
parser.add_argument('--outdir', type=str, required=False, help='Outdir for normalized file' )
parser.add_argument('--reference', type=str, required=True, help='Path to the Fasta file' )
parser.add_argument('--compress', type=bool, required=False, help='BGZIP compress the VCF output' )


args = parser.parse_args()

if __name__ == '__main__':
            
    vcfNorm = VcfNormalize(vcf=args.filename,vt_folder=args.vt_folder, bgzip_folder=args.bgzip_folder)
    vcfNorm.run_vtnormalize(outprefix=args.outprefix,reference=args.reference,compress=args.compress,outdir=args.outdir)
    
