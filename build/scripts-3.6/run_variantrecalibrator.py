from VcfFilter import VcfFilter
from ReseqTrackDB import File
from ReseqTrackDB import ReseqTrackDB

import argparse
import os
import logging
import datetime

#get command line arguments

parser = argparse.ArgumentParser(description='Script to run GATK VariantRecalibrator, which is part of the VQSR filtering procedure')

parser.add_argument('--gatk_folder', type=str, required=True, help='Folder containing GATK jar file' )
parser.add_argument('--vcf', type=str, required=True, help='Path to the VCF file that will be analysed' )
parser.add_argument('--reference', type=str, required=True, help='Path to the reference Fasta file' )
parser.add_argument('--resources', type=str, required=True, help='JSON File with resources that will be used by VariantRecalibrator' )
parser.add_argument('--outprefix', type=str, required=False, help='Prefix used for outputfiles' )
parser.add_argument('--intervals', type=str, required=False, help='i.e. chr20:10000-15000' )

args = parser.parse_args()

if __name__ == '__main__':

    VcfFilter = VcfFilter(vcf=args.vcf,caller='UG',gatk_folder=args.gatk_folder, reference=args.reference)

    outdict=VcfFilter.run_variantrecalibrator(args.resources,mode='SNP',outprefix=args.outprefix, intervals=args.intervals)



