with open("ggvp_samples_metadata.f2.txt") as f:
 for line in f:
  line=line.rstrip('\n')
  elms=line.split("\t")
  if elms[2]=="RELATED":
   print(elms[0]+"\tCHILD")
  else:
   print(elms[0]+"\t"+elms[1])
