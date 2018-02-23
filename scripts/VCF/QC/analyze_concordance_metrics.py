from ReseqTrackDB import *
from collections import defaultdict

import argparse
import pdb
import pandas as pd
import numpy as np
import pdb
import matplotlib.pyplot as plt
plt.switch_backend('agg') 
import numpy as np

parser = argparse.ArgumentParser(description='Script to perform the analysis of the results produced as a result of running\
 the PyHive::PipeConfig::CalcGtConcordance pipeline')

#RESEQTRACK DB conn params
parser.add_argument('--hostname', type=str, required=True, help='Hostname for ReseqTrack DB' )
parser.add_argument('--username', type=str, required=True, help='User for ReseqTrack DB' )
parser.add_argument('--port', type=int, required=True, help='Port number in the ReseqTrack DB' )
parser.add_argument('--pwd', type=str, help='PWD for the ReseqTrack DB' )
parser.add_argument('--db', type=str, required=True, help='DB name in the ReseqTrack DB' )

#Rest of params
parser.add_argument('--ifile', type=str, required=True, help='Path to VCF file that needs to be tracked in a ReseqTrack DB\
 used in the run of the PyHive::PipeConfig::CalcGtConcordance pipeline' )
parser.add_argument('--outprefix', type=str, required=True, help='Path to output file' )
parser.add_argument('--title', type=str, required=True, help='Title for plots' )
parser.add_argument('--pattern', type=str, required=True, help='If there is more than 1 set of attributes for the same file, then use the pattern in this arg to filter' )
parser.add_argument('--width', type=str, required=True, help='Width for plots. i.e. 600' )
parser.add_argument('--height', type=str, required=True, help='Height for plots. i.e. 25' )
parser.add_argument('--size', type=str, required=True, help='Size for font in plots. i.e. 50' )

args = parser.parse_args()

hostname=args.hostname
username=args.username
db=args.db
port=args.port
pwd=args.pwd

reseqtrackdb_object=ReseqTrackDB(host=hostname,user=username,port=port,pwd=pwd,db=db)

vcf_f=reseqtrackdb_object.fetch_file_by_url(args.ifile)

attrbs=reseqtrackdb_object.fetch_attributes_by_other_id(vcf_f.dbID)

def rec_dd():
    return defaultdict(rec_dd)


#valid attributes, edit it to modify the order of columns
valid_attrb=['VAR_SENSITIVITY','VAR_SPECIFICITY','GENOTYPE_CONCORDANCE']
x = rec_dd()

seen=0

for a in attrbs:
    astring=re.sub('^GT_CONC_','',a.name)
    astring=re.sub('LCvsAFFY_HASPED','LCvsAFFYHASPED',astring)
    astring=re.sub('LCvsAFFY_NOPED','LCvsAFFYNOPED',astring)
    astring=re.sub('EXvsAFFY_HASPED','EXvsAFFYHASPED',astring)
    astring=re.sub('LC_tsSNPS99_9_tsINDEL80_vsOMNI','LCtsSNPS999tsINDEL80vsOMNI',astring)
    astring=re.sub('LC_tsSNPS99_9_tsINDEL80_vsAFFY_HASPED','LCtsSNPS999tsINDEL80vsAFFYHASPED',astring)

    bits=astring.split("_")
    sample=bits.pop(0)
    prefix=bits.pop(0)
    if prefix==args.pattern:
        seen=1
        attrb="_".join(bits)
        if attrb in valid_attrb:
            x[sample][attrb][a.value]
if seen==0:
    print("No attrbs retrieved from DB using pattern {0}:".format(args.pattern))
    raise Exception
    
data=[]
for sample in x.keys():
    row=[]
    row.append(sample)
    for attr in valid_attrb:
        for value in x[sample][attr]:
            row.append(value)
    data.append(row)

df=pd.DataFrame(data,columns=['sample','VAR_SENSITIVITY','VAR_SPECIFICITY'
                              ,'GENOTYPE_CONCORDANCE'])
df.set_index('sample', inplace=True)

width=int(args.width)
height=int(args.height)
size=int(args.size)


#print stats
print("Median:")
for i,v in df.median().iteritems(): print(i,v)
print("\n")

def get_whiskers(a):
    upper_quartile = np.percentile(a, 75)
    lower_quartile = np.percentile(a, 25)
    iqr = upper_quartile - lower_quartile
    upper_whisker = a[a<=upper_quartile+1.5*iqr].max()
    lower_whisker = a[a>=lower_quartile-1.5*iqr].min()
    return (upper_whisker,lower_whisker)

(sens_uwhisker,sens_lwhisker)=get_whiskers(df['VAR_SENSITIVITY'].values)
(gconc_uwhisker,gconc_lwhisker)=get_whiskers(df['GENOTYPE_CONCORDANCE'].values)

print("Sensitivity uwhisker:{0} lwhisker:{1}".format(sens_uwhisker,sens_lwhisker))
print("Genotype concordance uwhisker:{0} lwhisker:{1}".format(gconc_uwhisker,gconc_lwhisker))

#plot sensitivity
plt.figure(1)
ax1=df['VAR_SENSITIVITY'].plot(rot=90,figsize=(width,height))
ax1.set_xticks(np.arange(len(df.index)))
ax1.set_xticklabels(df.index)

plt.ylabel("sensitivity", size=size)
plt.title(args.title, size=size)
plt.xlabel("samples", size=size)

fig1 = ax1.get_figure()
fig1.savefig('{0}_sensitivity.pdf'.format(args.outprefix),format='pdf')

#plot specificity
plt.figure(2)
ax2=df['VAR_SPECIFICITY'].plot(rot=90,figsize=(width,height))
ax2.set_xticks(np.arange(len(df.index)))
ax2.set_xticklabels(df.index)

plt.ylabel("specificity", size=size)
plt.title(args.title, size=size)
plt.xlabel("samples", size=size)

fig2 = ax2.get_figure()
fig2.savefig('{0}_specificity.pdf'.format(args.outprefix),format='pdf')

#plot genotype concordance
plt.figure(3)
ax3=df['GENOTYPE_CONCORDANCE'].plot(rot=90,figsize=(width,height))
ax3.set_xticks(np.arange(len(df.index)))
ax3.set_xticklabels(df.index)

plt.ylabel("genotype concordance", size=size)
plt.title(args.title, size=size)
plt.xlabel("samples", size=size)

fig3 = ax3.get_figure()
fig3.savefig('{0}_genotype_concordance.pdf'.format(args.outprefix),format='pdf')

#boxplot
df.plot.box(rot=90,figsize=(15,20),fontsize=10)
plt.title(args.title,size=20)
plt.savefig('{0}.boxplot.pdf'.format(args.outprefix),format='pdf')

