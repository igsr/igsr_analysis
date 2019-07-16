"""
This script will generate a manifest file (that can be used with nextflow download pipeline)
from an index file

usage: print get_manifest_from_index.py <aln.index> </path/to/download/>
"""

import sys
import os

print("url,dest")

with open(sys.argv[1]) as f:
 for line in f:
     line=line.rstrip("\n")
     if not line.startswith("#"):
         elms=line.split("\t")
         url=elms[0]
         basename=os.path.basename(url)
         dest="{0}/{1}".format(sys.argv[2],basename)
         print("{0},{1}".format(url,dest))
