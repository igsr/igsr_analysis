"""
This script will generate a manifest file (that can be used with nextflow download pipeline)
from an index file

usage: print get_manifest_from_index.py <aln.index>
"""
import sys
import os

print("url,name")

with open(sys.argv[1]) as f:
 for line in f:
     line=line.rstrip("\n")
     if not line.startswith("#"):
         elms=line.split("\t")
         url=elms[0]
         basename=os.path.basename(url)
         print("{0},{1}".format(url,basename))
