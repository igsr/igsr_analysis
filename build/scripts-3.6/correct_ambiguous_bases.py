#script to correct the ambiguous bases in an uncompressed vcf
import re
import sys

ref_count=0
alt_count=0

if (len(sys.argv)<2):
    sys.exit("[USAGE] python correct_ambiguous_bases.py in.vcf out.vcf")

infile=sys.argv[1]
outfile=sys.argv[2]

f=open(outfile,'w');

with open(infile,'r') as fin:
    for line in fin:
        if not line.startswith("#"):
            bits=line.split("\t")
            ref=bits[3]
            alt=bits[4]
            if re.search(r"[^ATGC.,]", ref):
                ref_count+=1
                ref=re.sub('[^ACGT.]','N',ref)
            if re.search(r"[^ATGC.,]", alt):
                alt_count+=1
                alt=re.sub('[^ACGT.]','N',alt)
            bits[3]=ref
            bits[4]=alt
            nline='\t'.join(bits)
            f.write(nline)
        else:
            f.write(line)
f.close()

print("Sites with ambiguous bases in the REF column is:{0}".format(ref_count))
print("Sites with ambiguous bases in the ALT column is:{0}".format(alt_count))
