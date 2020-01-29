with open("ggvp_samples_metadata.txt") as f:
 for line in f:
  line=line.rstrip('\n')
  if line.startswith("ID"):
   continue
  elms=line.split("\t")
  if elms[2]!="NA" or elms[2]!="NA":
   print(elms[0]+"\t"+elms[1]+"\tRELATED")
  else:
   print(elms[0]+"\t"+elms[1]+"\tUNRELATED")
