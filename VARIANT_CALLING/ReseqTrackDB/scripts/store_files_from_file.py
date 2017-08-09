'''

This script will parse a text file containing file metadata and will store each of the files described in the text file in a ReseqTrack DB

Created on 30 Nov 2016

@author: ernesto

'''

import argparse
from ReseqTrackDB import ReseqTrackDB
from ReseqTrackDB import File

#get command line arguments

#RESEQTRACK DB conn params
parser = argparse.ArgumentParser(description='This script will parse a text file containing file metadata and will store each of the files described in the text file in a ReseqTrack DB')
parser.add_argument('--hostname', type=str, required=True, help='Hostname for ReseqTrack DB' )
parser.add_argument('--username', type=str, required=True, help='User for ReseqTrack DB' )
parser.add_argument('--port', type=int, required=True, help='Port number in the ReseqTrack DB' )
parser.add_argument('--pwd', type=str, help='PWD for the ReseqTrack DB' )
parser.add_argument('--db', type=str, required=True, help='DB name in the ReseqTrack DB' )

parser.add_argument('--f', type=str, required=True, help='File with data on the files that want to be stored' )
'''
file passsed using the --f arg should have the following format:

<path>\t<md5><\t><type><\t><size><\t><host_id><\t><withdrawn><\t><created>

'''

args = parser.parse_args()

reseqdb = ReseqTrackDB(host=args.hostname,user=args.username,port=args.port,pwd=args.pwd,db=args.db)

with open(args.f) as f:
    for line in f:
        line=line.rstrip('\n')
        bits=line.split('\t')
        f=File(path=bits[0],type=bits[2],size=bits[3],md5=bits[1],host_id=bits[4],withdrawn=bits[5],created=bits[6])
        f.store(reseqdb)

