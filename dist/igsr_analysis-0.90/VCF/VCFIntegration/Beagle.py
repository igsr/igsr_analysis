'''
Created on 21 Jul 2017

@author: ernesto
'''
import os
import pdb
import subprocess
from Utils.RunProgram import RunProgram
from collections import namedtuple
import tempfile
import re

class Beagle(object):
    '''
    Class to operate on a VCF file and run Beagle and other Beagle-related operations on it
    '''

    def __init__(self, vcf, beagle_folder=None, beagle_jar=None, makeBGLCHUNKS_folder=None, prepareGenFromBeagle4_folder=None):
        '''
        Constructor

        Class variables
        ---------------
        vcf : str, Required
             Path to vcf file
        beagle_folder : str, Optional
                        Path to folder containing Beagle's jar file
        beagle_jar : str, Optional
                     Name of Beagle jar file. i.e. beagle.08Jun17.d8b.jar 
        makeBGLCHUNKS_folder : str, Optional
                               Path to folder containing makeBGLCHUNKS binary
                               (see https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gettingstarted)
        prepareGenFromBeagle4_folder : str, Optional
                                       Path to folder containing makeBGLCHUNKS binary
        '''

        if os.path.isfile(vcf) is False:
            raise Exception("File does not exist")

        self.vcf = vcf
        self.beagle_folder = beagle_folder
        self.beagle_jar = beagle_jar
        self.makeBGLCHUNKS_folder = makeBGLCHUNKS_folder
        self.prepareGenFromBeagle4_folder = prepareGenFromBeagle4_folder

    def run_beagle(self, outprefix, outdir=None, region=None, verbose=False, correct=False, **kwargs):
        '''
        Method that wraps Beagle (see https://faculty.washington.edu/browning/beagle/beagle.html)
        and will be used to call genotypes on a VCF file containing GT likelihoods

        Parameters
        ----------
        outprefix: str, required
              Prefix used for output file
        outdir : str, optional
                 outdir for output files
        region : str, optional
                 chr or chr interval that will be analyzed. i.e. chr20 or chr20:10000000-11000000
        verbose : bool, optional
                  if true, then print the command line used for running Beagle
        correct : bool, optional
                  Note: that it seems there is an incompatibility between zlib libraries used in Beagle4 and in BOOST on some platforms.
                  This involves either the last line of the file being skipped or a segfault. If correct=True, then this function will fix this issue
                  by recompressing the Beagle4 output files. Default=False
        window: int, optional
                number of markers to include in each sliding
                window. Default: 50000
        overlap: int, optional
                 specifies the number of markers of overlap between sliding
                 windows. Default: 3000
        niterations: unt, optional
                     specifies the number of phasing iterations. Default:
                     niterations=5
        nthreads : int, optional
                   number of threads. If not specified then the nthreads parameter 
                   will be set equal to the number of CPU cores on the host machine

        Returns
        -------
        Compressed VCF file with the genotype calls
        '''

        if self.beagle_folder is None or self.beagle_jar is None:
            raise Exception("Provide the folder containing the Beagle jar file and the Beagle jar file name")

        Arg = namedtuple('Argument', 'option value')
        args=[]

        outfile=""
        if outdir is not None:
            outfile="{0}/{1}.".format(outdir, outprefix)
        else:
            outfile="{0}.".format(outprefix)

        if region is not None:
            region_str=re.sub(":|-",".",region)
            outfile+="{0}.".format(region_str)
            args.append(Arg('chrom',region))

        outfile+="beagle"

        args.extend([Arg('gl',self.vcf), Arg('out',outfile)])

        for k,v in kwargs.items():
            args.append(Arg(k,v))

        runner=RunProgram(program='java -jar {0}/{1}'.format(self.beagle_folder,self.beagle_jar), args=args, arg_sep="=")

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout=runner.run_checkoutput()

        outfile=outfile+".vcf.gz"

        if correct is True:
            # creating temp file in order to perform the correction
            temp = tempfile.NamedTemporaryFile(delete=False)
            gzipRunner=RunProgram(program='gzip', parameters=['-c','>',temp.name])
            zcatRunner=RunProgram(program='zcat',parameters=[outfile],downpipe=[gzipRunner])

            if verbose is True:
                print("Command line for vcf.gz correction partA is: {0}".format(zcatRunner.cmd_line))
            
            #run zcat file | gzip -c > tmp.file 
            zcatRunner.run_checkoutput()

            #mv tmp.file back to outfile
            mvRunner=RunProgram(program='mv', parameters=[temp.name,outfile])

            if verbose is True:
                print("Command line for vcf.gz correction partB is: {0}".format(mvRunner.cmd_line))

            mvRunner.run_checkoutput()
            
        return outfile

    def make_beagle_chunks(self,window,overlap,outfile,verbose=False):
        '''
        Method to define chromosome chunks for Beagle
        see https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gettingstarted

        Parameters
        ----------
        window: int, required
                The chunk size (--window) in number of variant sites
        overlap: int, required
                 The overlap size (--overlap) in number of variant sites
        outfile: str, required
                 Name of output file. i.e. 'chunk.coordinates'
        verbose : bool, optional
                  if true, then print the command line used for running this tool.Default=False

        Returns
        -------
        Path to file with the coordinates of the chunk

        '''

        if self.makeBGLCHUNKS_folder is None:
            raise Exception("Provide the folder for the makeBGLCHUNKS binary")

        Arg = namedtuple('Argument', 'option value')

        args=[Arg('--vcf',self.vcf), Arg('--window',window), Arg('--overlap',overlap), Arg('--output',outfile)]

        runner=RunProgram(path="{0}/".format(self.makeBGLCHUNKS_folder), program='makeBGLCHUNKS', args=args)
        
        if verbose is True:
             print("Command line for running makeBGLCHUNKS is: {0}".format(runner.cmd_line))

        runner.run_checkoutput()
            
        return outfile

    def prepare_Gen_From_Beagle4(self,prefix_in,outprefix,threshold=0.995,verbose=False):
        '''
        Method that uses prepareGenFromBeagle4 in order to convert the different Beagle chunks
        generated by 'self.make_beagle_chunks' into a single concatenated output that can be used 
        with SHAPEIT.
        see https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gettingstarted

        Parameters
        ----------
        prefix_in: str, required
                   prefix used in the output of the different Beagle chunks after running  method 'self.run_beagle'.
                   i.e. output.beagle4.22.*.
        outprefix: str, required
                   Prefix used for output files. i.e. If prefix 'input.shapeit.chr22' is used. Then it will generate the following files:
                   input.shapeit.chr22.gen.gz
                   input.shapeit.chr22.gen.sample
                   input.shapeit.chr22.hap.gz
                   input.shapeit.chr22.hap.sample
        threshold: float, optional
                   Threshold meaning that all genotypes with a posterior above 0.995 are directly fixed and will only need phasing in the SHAPEIT step.
                   Default: 0.995
        verbose : bool, optional
                  if true, then print the command line used for running this tool.Default=False

        Returns
        -------
        A dict with the path to the 4 output files (*.gen.* and *.hap.*) that can be used with SHAPEIT
        
        '''

        if self.prepareGenFromBeagle4_folder is None:
            raise Exception("Provide the folder for the prepareGenFromBeagle4 binary")

        posteriors="{0}*.vcf.gz".format(prefix_in)

        Arg = namedtuple('Argument', 'option value')

        args=[Arg('--likelihoods',self.vcf), Arg('--posteriors',posteriors), Arg('--output',outprefix)]

        runner=RunProgram(path="{0}/".format(self.prepareGenFromBeagle4_folder), program='prepareGenFromBeagle4', args=args)

        if verbose is True:
             print("Command line for running prepareGenFromBeagle4 is: {0}".format(runner.cmd_line))

        runner.run_checkoutput()

        outdict={ 'gen_gz' :'{0}.gen.gz'.format(outprefix),
                  'gen_sample' : '{0}.gen.sample'.format(outprefix),
                  'hap_gz' : '{0}.hap.gz'.format(outprefix),
                  'hap_sample' : '{0}.hap.sample'.format(outprefix)
        }
        
        return outdict
