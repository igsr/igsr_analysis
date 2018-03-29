'''
Created on 24 Apr 2017

@author: ernesto
'''
import pdb
import re
import subprocess
import tempfile

class BEDTools(object):
    '''
    Class to represent different operations performed with BEDTools
    '''

    def __init__(self, bedtools_folder):
        '''
        Constructor

        Class variables
        ---------------
        bedtools_folder : str, Optional
                          Path to folder with bedtools binary
        '''

        self.bedtools_folder = bedtools_folder

    def make_windows(self, w, g, s=None, subtract=None, lextend=None, rextend=None, verbose=False):
        '''
        This method will make windows from a genome file by using 'bedtools makewindows'

        Parameters
        ----------
        w : int , Required
            width of windows in bp
        g : str, Required
            Path to genome file
        s : int , Optional
           overlap in bp. i.e. if -w 100 -s 80 will generate:

           chr1    0       100
           chr1    80      180
           chr1    160     260
           ...
           So, -s defines the offset in bp

           Another example -w 1000 -s 200

           chr1    0       1000
           chr1    200     1200
           chr1    400     1400
           chr1    600     1600
        lextend : int, Optional
                  Extend each interval to the left by int bases 

        rextend : int, Optional
                  Extend each interval to the right by int bases
        
        subtract : str, Optional
                   BED file containing the features that will be removed from the generated windows.
                   For example, if we have the following window:
                   
                   chr20 1000 2000

                   And we have the following feature in the BED file: chr20 1100 1200
                   Then the resulting windows will be like:
                     
                   chr20 1000 1100
                   chr20 1200 2000
        
        verbose : boolean, optional.
                  Default=False

        Returns
        ------
        A list of lists. Each sublist is composed of ['chr','start','end']
        
        It will return an empty list if not elements for a certain chr are defined
        '''

        command = ""
        if self.bedtools_folder:
            command += self.bedtools_folder+"/"

        command += "bedtools makewindows -g {0} -w {1}".format(g,w)

        if s is not None:
            command += " -s {0}".format(s)

        coordlist=[]

        if verbose is not False:
            print(command)

        try:
            stdout=subprocess.check_output(command, shell=True)
            coordlist=[l.split("\t") for l in stdout.decode("utf-8").strip().split("\n")]
        except subprocess.CalledProcessError as exc:
            raise Exception(exc.output)

        if subtract is not None:
            temp = tempfile.NamedTemporaryFile()
            try:
                f=open(temp.name,'w');
                for i in coordlist:
                    f.write("{0}\t{1}\t{2}\n".format(i[0],i[1],i[2]))
                f.close()
                
                command1= "{0}/bedtools subtract -a {1} -b {2}".format(self.bedtools_folder, 
                                                                       temp.name, subtract)
                coordlist=None
                try:
                    stdout=subprocess.check_output(command1, shell=True)
                    coordlist=[l.split("\t") for l in stdout.decode("utf-8").strip().split("\n")]
                except subprocess.CalledProcessError as exc:
                    raise Exception(exc.output)
            finally:
                temp.close()

        if lextend is not None:
            first_seen=False
            for k,l in enumerate(coordlist):
                if first_seen==True: l[1]=str(int(l[1])+lextend)
                first_seen=True
                coordlist[k]=l

        if rextend is not None:
            for k,l in enumerate(coordlist):
                if k!=len(coordlist)-1: l[2]=str(int(l[2])+rextend)
                coordlist[k]=l

        return coordlist

    def __str__(self):
        sb = []
        for key in self.__dict__:
            sb.append("{key}='{value}'".format(key=key, value=self.__dict__[key]))

        return ', '.join(sb)

    def __repr__(self):
        return self.__str__()
