'''
Created on 31 Mar 2017

@author: ernesto
'''

import os
import subprocess

class VcfNormalize(object):
    '''
    Normalize variants in a VCF file
    '''


    def __init__(self, vcf, vt_folder=None, vcflib_folder=None, bgzip_folder=None,
                 gatk_folder=None):
        '''
        Constructor

        Class variables
        ---------------
        vcf : str, Required
             Path to gzipped vcf file
        vt_folder : str, Optional
                      Path to folder containing the vt binary
        vcflib_folder : str, Optional
                        Path to folder containing the different vcflib binaries
        bgzip_folder : str, Optional
                       Path to folder containing the bgzip binary
        gatk_folder : str, Optional
                      Path to folder containing the GATK jar file
        '''

        if os.path.isfile(vcf) is False:
            raise Exception("File does not exist")

        self.vcf = vcf
        self.vt_folder = vt_folder
        self.vcflib_folder = vcflib_folder
        self.bgzip_folder = bgzip_folder
        self.gatk_folder = gatk_folder

    def run_vtnormalize(self, outprefix, reference, compress=False, outdir=None):
        '''
        Run vt normalize on a vcf file

        Parameters
        ----------
        outprefix : str, required
              prefix for outputfile
        reference : str, required
              path to Fasta file with reference
        compress : boolean, optional
              bgzip compress the normalized VCF
        outdir : str, optional
            If provided, then put output files in this folder

        Returns
        -------
        A string with path to normalized file
        '''

        command = ""
        if self.vt_folder:
            command += self.vt_folder+"/"

        if outdir:
            outprefix = "{0}/{1}".format(outdir, outprefix)

        outprefix = outprefix+".norm.vcf"

        if compress is True:
            outprefix = outprefix+".gz"
            bgzip_path = ""
            if self.bgzip_folder:
                bgzip_path = "{0}/bgzip".format(self.bgzip_folder)
            else:
                bgzip_path = "bgzip"

            command += "vt normalize -m -r {0} {1}  | {2} -c > {3} ".\
            format(reference, self.vcf, bgzip_path, outprefix)

        elif compress is None or compress is False:
            command += "vt normalize -m -r {0} -o {1} {2}".format(reference, outprefix, self.vcf)

        try:
            subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)
        except subprocess.CalledProcessError as exc:
            print(exc.output)
            print('ERROR RUNNING COMMAND: {0} '.format(exc.cmd))

        return outprefix

    def run_vcfallelicprimitives(self, outprefix, compress=False, outdir=None,
                                 keepinfo=True, keepgeno=True):
        '''
        Run vcfallelicprimitives on a vcf file

        This program is used to decompose complex variants into a canonical SNP and
        indel representation,generating phased genotypes for available samples.

        Parameters
        ----------

        outprefix : str, required
              prefix for outputfiles
        compress : boolean, optional
              bgzip compress the normalized VCF
        outdir : str, optional
            If provided, then put output files in this folder
        keepinfo : bool, optional. Default=True
            Maintain site and allele-level annotations when decomposing.
            Note that in many cases, such as multisample VCFs, these won't
            be valid post-decomposition.  For biallelic loci in single-sample
            VCFs, they should be usable with caution
        keepgeno : bool, optional. Default=True
            Maintain genotype-level annotations when decomposing.  Similar
            caution should be used for this as for keep-info.

        Returns
        -------
        A string with path to decomposed file
        '''

        command = ""
        if self.vcflib_folder:
            command += self.vcflib_folder+"/"

        if outdir: 
            outprefix = "{0}/{1}".format(outdir, outprefix)

        outprefix = outprefix+".aprimitives.vcf"

        command += "vcfallelicprimitives {0} ".format(self.vcf)

        if keepinfo is True:
            command += "--keep-info "

        if keepgeno is True:
            command += "--keep-geno "

        if compress is True:
            outprefix = outprefix+".gz"
            bgzip_path = ""
            if self.bgzip_folder:
                bgzip_path = "{0}/bgzip".format(self.bgzip_folder)
            else:
                bgzip_path = "bgzip"

            command += "| {0} -c > {1} ".format(bgzip_path, outprefix)

        elif compress is None or compress is False:
            command += " > {0}".format(outprefix)

        try:
            subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)
        except subprocess.CalledProcessError as exc:
            print(exc.output)
            print('ERROR RUNNING COMMAND: {0} '.format(exc.cmd))

        return outprefix

    def run_gatk_VariantsToAllelicPrimitives(self, outprefix, reference,
                                             outdir=None, compress=None):
        '''
        Run GATK VariantsToAllelicPrimitives in order to decompose MNPs
         into more basic/primitive alleles

        Parameters
        ----------

        outprefix : str, required
                    prefix for outputfiles
        reference : str, Required
                     Path to fasta file containing the reference
        outdir : str, optional
                   If provided, then put output files in this folder
        compress : boolean, optional
                   bgzip compress the normalized VCF

        Returns
        -------
        A string with path to decomposed file

        '''

        if self.gatk_folder is None:
            raise Exception("Error. I need that the folder containing the GATK "
                            "jar file is defined!")

        if outdir: 
            outprefix = "{0}/{1}".format(outdir, outprefix)

        outprefix = outprefix+".aprimitives.vcf"

        command = "java -jar {0}/GenomeAnalysisTK.jar -T VariantsToAllelicPrimitives -R {1} -V {2}"\
        " -o {3}".format(self.gatk_folder, reference, self.vcf, outprefix)

        try:
            subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            print(exc.output)
            print('ERROR RUNNING COMMAND: {0} '.format(exc.cmd))

        if compress == 'True':

            bgzip_path = ""
            if self.bgzip_folder:
                bgzip_path = "{0}/bgzip".format(self.bgzip_folder)
            else:
                bgzip_path = "bgzip"

            command = "{0} -c {1} > {1}.gz ".format(bgzip_path, outprefix)

            try:
                subprocess.check_output(command, shell=True)
            except subprocess.CalledProcessError as exc:
                print(exc.output)
                print('ERROR RUNNING COMMAND: {0} '.format(exc.cmd))

            #delete tmp files
            os.remove(outprefix)
            os.remove(outprefix+".idx")
            outprefix += ".gz"
            return outprefix

        elif compress == 'False' or compress is None:
            return outprefix
        else:
            raise Exception("Compress option: {0} is not recognized".format(compress))
