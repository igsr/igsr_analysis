'''
Created on 03 Aug 2017

@author: ernesto
'''

import os
import subprocess
import pdb
import gzip
import re
from Utils.RunProgram import RunProgram
from collections import namedtuple

class VcfUtils(object):
    '''
    Class to represent a misc of actions that can be done on a single or multiple VCF files
    '''


    def __init__(self, vcf=None, vcflist=None, bcftools_folder=None, bgzip_folder=None, gatk_folder=None, java_folder=None):
        '''
        Constructor

        Class variables
        ---------------
        vcf : str, optional
             Path to gzipped vcf file
        vcflist : list, optional
             List of dicts containing setname:vcf_paths (keys:values) pairs
        bcftools_folder : str, Optional
                         Path to folder containing the bcftools binary
        bgzip_folder : str, Optional
                       Path to folder containing the bgzip binary
        gatk_folder : str, Optional
                      Path to folder containing the jar file
        java_folder : str, Optional
                      Path to folder containing the java binary

        Imp: Either 'vcf' or 'vcflist' variables should be initialized
        '''

        if not vcf and not vcflist:
            raise Exception("Either a vcf file or a list of vcf files should be used\
                             to initialize this class")

        if vcf is not None:
            if os.path.isfile(vcf) is False:
                raise Exception("File does not exist")

        self.vcf = vcf
        self.vcflist = vcflist
        self.bcftools_folder = bcftools_folder
        self.bgzip_folder = bgzip_folder
        self.gatk_folder = gatk_folder
        self.java_folder = java_folder

    def reheader(self, newheader, outprefix, samplefile=None, verbose=False):
        '''
        Modifiy the VCF's header with the newheader

        Parameters
        ----------
        newheader : string, required
                   Path to the file containing the new header
        outprefix : string, required
                    prefix for output files
        samplefile : string, optional
                     Path to the file with the sample names that will included
                     in the new header
        verbose : bool, optional
                  increase the verbosity, default=False

        Returns
        -------
        Path to the VCF with the modified header
        '''
    
        outfile=outprefix+".reheaded.vcf.gz"

        Arg = namedtuple('Argument', 'option value')

        args=[Arg('-h', newheader), Arg('-o',outfile)]

        if samplefile is not None:
            args.append(Arg('-s',samplefile))

        runner=RunProgram(path=self.bcftools_folder, program='bcftools reheader', args=args, parameters=[self.vcf])

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout=runner.run_checkoutput()

        return outfile

    def combine(self, labels, reference, outprefix, compress=False, outdir=None, ginterval=None, 
                genotypemergeoption=None, filteredrecordsmergetype=None, threads=1, options=None, verbose=False):
        '''
        Combine VCFs using GATK's CombineVariants into a single VCF

        Parameters
        ----------
        labels : list, required
                 List of labels used for each of the VCFs in self.vcflist. The order of the labels
                 should be the same that the VCFs in the list
        reference : str, required
                    Path to Fasta file with reference
        outprefix : str, required
                    prefix used for output file
        compress : boolean, optional
                   Compress the output VCF with bgzip. Default=False
        outdir : str, optional
                 Path to folder used to write the results to
        ginterval : str, optional
                    Genomic interval used to restrict the analysis. i.e. chr20:1000-2000
        genotypemergeoption : str, optional
                    Determines how we should merge genotype records for samples shared across the ROD files. 
                    Possible values are: UNIQUIFY, PRIORITIZE, UNSORTED, REQUIRE_UNIQUE
        filteredrecordsmergetype : str, optional
                    Determines how we should handle records seen at the same site in the VCF, but with different FILTER fields
                    Possible values are : KEEP_IF_ANY_UNFILTERED, KEEP_IF_ANY_UNFILTERED, KEEP_UNCONDITIONAL
        threads : int, optional
                  Number of trades to use. Default=1
        options : list, optional
                   List of options. i.e. ['-env','--filteredAreUncalled']
        verbose : bool, optional
                  increase the verbosity, default=False
    
        Returns
        -------
        Path to the merged VCF
        '''
        
        Arg = namedtuple('Argument', 'option value')

        args=[Arg('-T','CombineVariants'), Arg('-R',reference), Arg('-nt', threads)]

        variants_str=""
        for path, label in zip(self.vcflist, labels):
            if os.path.isfile(path) == False:
                print("Error reading from {0}".format(path))
                raise Exception("File does not exist")
            args.append(Arg('-V:{0}'.format(label),path))

        outfile=""
        if outdir:
            outfile= "{0}/".format(outdir)
        outfile+= "{0}.vcf".format(outprefix)
        
        if ginterval is not None:
            args.append(Arg('-L',ginterval))

        if genotypemergeoption is not None:
            args.append(Arg('--genotypemergeoption', genotypemergeoption))

        if filteredrecordsmergetype is not None:
            args.append(Arg('--filteredrecordsmergetype', filteredrecordsmergetype))

        params=[]
        if options:
            for opt in options:
                params.append(opt)

        pipelist=None
        if compress is True:
            outfile += ".gz"
            compressRunner=RunProgram(path=self.bgzip_folder,program='bgzip',parameters=[ '-c', '>', outfile])
            pipelist=[compressRunner]
        else:
            args.append(Arg('-o', outfile))

        runner=RunProgram(path=self.java_folder, program='java -jar {0}/GenomeAnalysisTK.jar'.format(self.gatk_folder), args=args, parameters=params, downpipe=pipelist)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout,stderr=runner.run_popen()

        return outfile


    def rename_chros(self, chr_types, outfile, compress=True):
        '''
        Function to modify the chr names in the VCF file
        For example:
        If file has UCSC-type chr names (i.e. chr1,chr2,...) then this 
        function will convert the UCSC-type chr names to Ensembl-type 
        chr names (i.e. 1,2,...) or vice-versa

        Parameters
        ----------
        chr_types : string, required
                    Type of chr names that will be written to the file. 
                    Possible values are: 
                         'ucsc'/'ensembl'
        outfile : string, required
                  File used for the output VCF
        compress : boolean, optional
                   Default: True
        
        Returns
        -------
        Path to the VCF with the chrosomes renamed
        
        '''

        command=""
        if chr_types=='ensembl':
            if compress is True:
                command += "zcat {0} | awk '{{gsub(/^chr/,\"\"); print}}' - | {1}/bgzip -c > {2}".format(self.vcf,self.bgzip_folder, outfile)
            else:
                command += "zcat {0} | awk '{{gsub(/^chr/,\"\"); print}}' - > {1}".format(self.vcf, outfile)
        elif chr_types=='ucsc':
            if compress is True:
                command += "zcat {0} | awk '{{if($0 !~ /^#/) print \"chr\"$0; else print $0}}' - | {1}/bgzip -c > {2}".format(self.vcf,self.bgzip_folder, outfile)
            else:
                command += "zcat {0} | awk '{{if($0 !~ /^#/) print \"chr\"$0; else print $0}}' - > {1}".format(self.vcf, outfile)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            print("Something went wrong.\n"
                  "Command used was: %s" % command)
            raise Exception(exc.output)



        return outfile

    def correct_ambiguity_codes(self,outfile):
        '''
        Function to correct the ambiguity bases in the VCF. This ambiguity
        may appear in the REF or ALT columns

        Parameters
        ----------
        outfile : string, required
                  File where the output VCF will be written

        Returns
        -------
        Path to vcf.gz file compressed with GZIP
        '''
        
        ref_count=0
        alt_count=0

        f=gzip.open(outfile,'wb');

        with gzip.open(self.vcf,'r') as fin:
            for line in fin:
                if not line.startswith(b"#"):
                    bits=line.split(b"\t")
                    ref=bits[3].decode("utf-8")
                    alt=bits[4].decode("utf-8")
                    if re.search(r"[^ATGC.,]", ref):
                        ref_count+=1
                        ref=re.sub('[^ACGT.]','N',ref)
                    if re.search(r"[^ATGC.,]", alt):
                        alt_count+=1
                        alt=re.sub('[^ACGT.]','N',alt)
                    bits[3]=ref.encode('utf-8')
                    bits[4]=alt.encode('utf-8')
                    nline=b'\t'.join(bits)
                    f.write(nline)
                else:
                    f.write(line)
        f.close()
        
        print("Sites with ambiguous bases in the REF column is:{0}".format(ref_count))
        print("Sites with ambiguous bases in the ALT column is:{0}".format(alt_count))

        return outfile

    def drop_genotypes(self, outfile, verbose=False):
        '''
        Function to drop the Genotype information from a VCF.
        This function uses bcftools -G to perform this operation

        Parameters
        ----------
        outfile : string, required
                  File where the output VCF will be written
        verbose : bool, optional
                  increase the verbosity, default=False

        Returns
        -------
        Path to the vcf.gz file without the GT information
        '''
        
        Arg = namedtuple('Argument', 'option value')
         
        args=[Arg('-o',outfile), Arg('-O', 'z')]

        runner=RunProgram(path=self.bcftools_folder, 
                          program='bcftools view -G', args=args, parameters=[self.vcf])

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout=runner.run_checkoutput()

        return outfile

    def drop_info(self, outfile, verbose=False):
        '''
        Function to remove the INFO annotation from a VCF.
        This function uses bcftools annotate  to perform this operation

        Parameters
        ----------
        outfile : string, required
                  File where the output VCF will be written
        verbose : bool, optional
                  increase the verbosity, default=False

        Returns
        -------
        Path to the vcf.gz file without the INFO annotation
        '''
        
        Arg = namedtuple('Argument', 'option value')

        args=[Arg('-o',outfile), Arg('-O','z')]

        runner=RunProgram(path=self.bcftools_folder,
                          program='bcftools annotate --remove INFO', args=args, parameters=[self.vcf])

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout=runner.run_checkoutput()

        return outfile

    def convert_PL2GL(self, outfile, threads=1, verbose=False):
        '''
        Function to convert PL fields into GT.
        This function makes use of Bcftools +tag2tag plugin
        
        Parameters
        ----------
        outfile : string, required
                  File where the output VCF will be written
                   Possible values are : KEEP_IF_ANY_UNFILTERED, KEEP_IF_ANY_UNFILTERED, KEEP_UNCONDITIONAL
        threads : int, optional
                  Number of trades to use. Default=1
        verbose : bool, optional
                  increase the verbosity, default=False

        Returns
        -------
        Path to the vcf.gz file with the PL fields converted
        
        '''

        Arg = namedtuple('Argument', 'option value')

        params=[self.vcf,'-Oz','--','-r','--pl-to-gl']

        runner=RunProgram(path=self.bcftools_folder,
                          program='bcftools +tag2tag', 
                          args=[Arg('--threads',threads), Arg('-o',outfile)], 
                          parameters=params)

        if verbose is True:
            print("Command line is: {0}".format(runner.cmd_line))

        stdout=runner.run_checkoutput()

        return outfile
        
        
