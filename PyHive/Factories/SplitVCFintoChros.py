import eHive
import os
import sys
from VCFfilter.BCFTools import BCFTools

class SplitVCFintoChros(eHive.BaseRunnable):
    """Split a VCF into the chromosomes present in a Fasta index"""

    def param_defaults(self):
        return {
        }

    def __str_to_bool(self,s):
        if s == 'True':
            return True
        elif s == 'False':
            return False
        else:
            raise ValueError # evil ValueError!

    def run(self):

        self.warning('Split file: %s'% self.param_required('filepath'))
        
        ifile=self.param_required('filepath')
        file=os.path.split(ifile)[1]
        outdir=self.param_required('work_dir')
        faix = self.param_required('faix')

        verbose=None
        if self.param_is_defined('verbose'):
            verbose=self.__str_to_bool(self.param('verbose'))

        bcftools_o = BCFTools(vcf=ifile,bcftools_folder=self.param('bcftools_folder'))

        files=[]
        ix=1
        for line in open(faix):
            if line.startswith("\n"):
                continue
            chr=line.split('\t')[0]
            self.warning('Splitting %s'% chr)
            chr_folder=outdir+"/"+chr
            if not os.path.isdir(chr_folder):
                os.mkdir(chr_folder)
            vcffile=bcftools_o.subset_vcf(region=chr,outprefix=file,outdir=chr_folder,
                                          create_index=True, threads=self.param('threads'), action='include', verbose=verbose)
            files.append(
                {
                    'chr': vcffile,
                    'ix': ix
                }
            )
            ix+=1

        self.param('files', files)

    def write_output(self):
        self.warning('{0} files have been created'.format(len(self.param('files'))))

        if self.param('verbose')=="True":
            for f in self.param('files'):
                self.warning("Chr file is %s" % f)
                
        self.dataflow(self.param('files'), 2)
