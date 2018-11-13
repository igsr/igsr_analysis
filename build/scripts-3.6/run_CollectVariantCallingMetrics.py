from VcfQC import VcfQC

import argparse
import os

#get command line arguments

parser = argparse.ArgumentParser(description='Script to run Picard CollectVariantCallingMetrics on a VCF file')

parser.add_argument('--picard_folder', type=str, required=True, help='Folder containing Picard jar file' )
parser.add_argument('--vcf', type=str, required=True, help='Path to the VCF file that will be analysed' )
parser.add_argument('--outprefix', type=str, required=True, help='Prefix for output files' )
parser.add_argument('--truth_vcf', type=str, required=True, help='VCF that will be used as truth callset' )
parser.add_argument('--interval_file', type=str, required=False, help='File with target intervals to restrict analysis to')

args = parser.parse_args()

if __name__ == '__main__':
    VcfQC = VcfQC(vcf=args.vcf,picard_folder=args.picard_folder)
    cvcmetrics=VcfQC.run_CollectVariantCallingMetrics(outprefix=args.outprefix,truth_vcf=args.truth_vcf,intervals=args.interval_file)

