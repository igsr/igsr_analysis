'''
Created on 21 Jul 2017

@author: ernesto
'''
import os
import pdb
import subprocess
import tempfile
import re

class Beagle(object):
    '''
    Class to operate on a VCF file and run Beagle and other Beagle-related operations on it
    '''

    def __init__(self, vcf, beagle_folder=None, makeBGLCHUNKS_folder=None, prepareGenFromBeagle4_folder=None):
        '''
        Constructor

        Class variables
        ---------------
        vcf : str, Required
             Path to vcf file
        beagle_folder : str, Optional
                        Path to folder containing Beagle's jar file
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

        program_folder = ""
        if self.beagle_folder:
            program_folder += self.beagle_folder + "/"

        outfile=""
        if outdir is not None:
            outfile="{0}/{1}.".format(outdir, outprefix)
        else:
            outfile="{0}.".format(outprefix)

        if region is not None:
            region_str=re.sub(":|-",".",region)
            outfile+="{0}.".format(region_str)

        outfile+="beagle"

        command = "java -jar {0}/beagle.08Jun17.d8b.jar gl={1} out={2}".format(program_folder,
                                                                               self.vcf,
                                                                               outfile)

        if region is not None:
            command += " chrom={0}".format(region)

        for k,v in kwargs.items():
            command += " {0}={1}".format(k,v)
            
        if verbose==True:
            print("Command used was: %s" % command)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            print("Something went wrong while running Beagle\n"
                  "Command used was: %s" % command)
            raise Exception(exc.output)
            
        outfile=outfile+".vcf.gz"

        if correct is True:
            # creating temp file in order to perform the correction
            temp = tempfile.NamedTemporaryFile(delete=False)
            try:
                correct_cmd1= "zcat {0} |gzip -c > {1} ".format(outfile,temp.name)
                subprocess.check_output(correct_cmd1, shell=True)
                correct_cmd2= "mv {0} {1}".format(temp.name,outfile)
                subprocess.check_output(correct_cmd2, shell=True)
            except subprocess.CalledProcessError as exc:
                print("Something went wrong while performing the segfault correction\n")
                raise Exception(exc.output)
            finally:
                # Automatically cleans up the file
                temp.close()

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

        program_folder = ""
        if self.makeBGLCHUNKS_folder:
            program_folder += self.makeBGLCHUNKS_folder + "/"

        command = "{0}/makeBGLCHUNKS --vcf {1} --window {2} --overlap {3} --output {4}".format(program_folder,
                                                                                                self.vcf,
                                                                                                window,
                                                                                                overlap,
                                                                                                outfile)
        
        if verbose==True:
            print("Command used was: %s" % command)

        try:
            subprocess.check_output(command, shell=True)
            if os.path.isfile(outfile) == False:
                print("Error. Something went wrong while running command: {0}".format(command))
                raise Exception("File cound not be created")
        except subprocess.CalledProcessError as exc:
            print("Something went wrong while running makeBGLCHUNKS\n"
                  "Command used was: %s" % command)
            raise Exception(exc.output)
        
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

        program_folder = ""
        if self.prepareGenFromBeagle4_folder:
            program_folder += self.prepareGenFromBeagle4_folder + "/"

        posteriors="{0}*.vcf.gz".format(prefix_in)

        command = "{0}/prepareGenFromBeagle4 --likelihoods {1} --posteriors {2} --threshold {3} --output {4}".format(program_folder,
                                                                                                                     self.vcf,
                                                                                                                     posteriors,
                                                                                                                     threshold,
                                                                                                                     outprefix)
    
        if verbose==True:
            print("Command used was: %s" % command)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            print("Something went wrong while running prepareGenFromBeagle4\n"
                  "Command used was: %s" % command)
            raise Exception(exc.output)

        outdict={ 'gen_gz' :'{0}.gen.gz'.format(outprefix),
                  'gen_sample' : '{0}.gen.sample'.format(outprefix),
                  'hap_gz' : '{0}.hap.gz'.format(outprefix),
                  'hap_sample' : '{0}.hap.sample'.format(outprefix)
        }
        
        return outdict
