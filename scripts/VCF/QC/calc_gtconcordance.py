'''
author: ernesto lowy (ernestolowy@gmail.com)

This script will parse a .tsv file in the following format:
POS     REF     ALT     GT
18948096        G       T       0|1
18948246        G       A       0|1
18948886        T       A       0|1
18949107        A       G       0|1
...

And will generate some contigency tables in order to have a feeling
on the Genotype concordance of the 2 call sets

USAGE: python calc_gtconcordance <igsr.tsv> <giab.tsv> 
'''
import pandas as pd
import sys

#read-in the igsr tsv file
DF_igsr=pd.read_csv(sys.argv[1], sep='\t',index_col=0)

#read-in the GIAB tsv file
DF_giab=pd.read_csv(sys.argv[2], sep='\t',index_col=0)

#combine 2 dataframes by index
final_DF=pd.merge(DF_igsr,DF_giab,left_index=True, right_index=True, suffixes=('_igsr', '_giab'))

# contigency tables for ref and alt alleles
ref_table=pd.crosstab(final_DF.REF_igsr, final_DF.REF_giab)
print("\n###Correspondence between REF alleles:")
print(ref_table)

alt_table=pd.crosstab(final_DF.ALT_igsr, final_DF.ALT_giab)
print("\n###Correspondence between ALT alleles:")
print(alt_table)

#contigency tables for genotype concordance
gt_tables=pd.crosstab(final_DF.GT_igsr, final_DF.GT_giab,margins=True, margins_name="Total")
print("\n### GT concordance (counts):")
print(gt_tables)

gt_tables=pd.crosstab(final_DF.GT_igsr, final_DF.GT_giab,margins=True, margins_name="Total",normalize=True).round(4)*100
print("\n### GT concordance (%):")
print(gt_tables)

