from Coord import *
import os

if __name__ == '__main__':
    
    fai_file="/Users/ernesto/Documents/workspace/Coordinates/files/only3chrs.fa.fai"
    
    if os.path.isfile(fai_file) == False:
        raise Exception("File does not exist")
    
    with open(fai_file) as f:
        for line in f: 
            elms=line.split("\t")
            c=Coord(id=elms[0],start=0,end=int(elms[1]))
            l=c.make_windows(step=10000000)
