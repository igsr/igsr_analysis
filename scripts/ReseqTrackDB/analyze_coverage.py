'''

This script will create a Boxplot with the median coverage calculated over a certain BAM file

Created on 01 Dec 2016

@author: ernesto

'''

import pandas as pd
import argparse
import matplotlib
matplotlib.use('Agg')
from ReseqTrackDB import ReseqTrackD

#get command line arguments

#RESEQTRACK DB conn params
parser = argparse.ArgumentParser(description='This script will create a Boxplot with the median coverage calculated over a certain BAM file')
parser.add_argument('--hostname', type=str, required=True, help='Hostname for ReseqTrack DB' )
parser.add_argument('--username', type=str, required=True, help='User for ReseqTrack DB' )
parser.add_argument('--port', type=int, required=True, help='Port number in the ReseqTrack DB' )
parser.add_argument('--pwd', type=str, help='PWD for the ReseqTrack DB' )
parser.add_argument('--db', type=str, required=True, help='DB name in the ReseqTrack DB' )

parser.add_argument('--type', type=str, required=True, help='Type in the file table that will be analysed' )
parser.add_argument('--outdir', type=str, required=True, help='Folder to write the pdf file with boxplots' )

args = parser.parse_args()

reseqdb = ReseqTrackDB(host=args.hostname,user=args.username,port=args.port,pwd=args.pwd,db=args.db)

files=reseqdb.fetch_files_by_type(args.type)

medians=[]
means=[]

for f in files:
    attbs=reseqdb.fetch_attributes_by_other_id(f.dbID)
    means.append(list(filter(lambda x:'PICARD:MEDIAN_COVERAGE' in x.name, attbs))[0].value)
    medians.append(list(filter(lambda x:'PICARD:MEAN_COVERAGE' in x.name, attbs))[0].value)


data=list(zip(means,medians))

df = pd.DataFrame(data = data, columns=['means','medians'])

#means
print("Mean of coverages\n=================\n")
print(df['means'].describe())

#medians
print("\nMedian of coverages\n===================\n")
print(df['medians'].describe())

props = dict(boxes="lightGreen", whiskers="DarkOrange", medians="DarkBlue", caps="Gray")

ax=df.plot.box(grid=True,return_type='axes',color=props, patch_artist=True, title="Coverage")


fig = ax.get_figure()
fig.savefig(args.outdir+"/BAM_coverage_boxplots.pdf",format='pdf')
