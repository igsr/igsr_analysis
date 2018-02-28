'''
Created on 24 Apr 2017

@author: ernesto
'''
import pdb
import re
import subprocess

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

    def make_windows(self, w, g, s=None, verbose=False):
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

        return coordlist

    def __str__(self):
        sb = []
        for key in self.__dict__:
            sb.append("{key}='{value}'".format(key=key, value=self.__dict__[key]))

        return ', '.join(sb)

    def __repr__(self):
        return self.__str__()
