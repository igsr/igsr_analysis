import vcf
import re
import sys

if (len(sys.argv)<2):
    sys.exit("[USAGE] python correct_ambiguous_bases.py in.vcf.gz out.vcf")

infile=sys.argv[1]
outfile=sys.argv[2]

vcf_reader = vcf.Reader(filename=infile)
vcf_writer = vcf.Writer(open(outfile, 'w'), vcf_reader)

#number of ambiguous bases in the REF/ALT column
no_ref=0
no_alt=0

for record in vcf_reader:
    ref=record.REF
    if re.search(r"[^ATGC.,]", ref):
        ref=re.sub('[^ACGT.]','N',ref)
        no_ref=no_ref+ref.count("N")

    alt=record.ALT
    new_alt=[]
    for i in alt:
        i=str(i)
        if re.search(r"[^ATGC.,]", i):
            i=re.sub('[^ACGT.]','N',i)
            no_alt=no_alt+i.count("N")
        new_alt.append(i)
    record.REF=ref
    record.ALT=new_alt
    vcf_writer.write_record(record)

print("Number of ambiguous bases in the REF column:{0}".format(no_ref))
print("Number of ambiguous bases in the ALT column:{0}".format(no_alt))
