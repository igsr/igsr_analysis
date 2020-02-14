import eHive
import os
from VCF.VCFfilter.BCFTools import BCFTools

class SubsetVcfWithBed(eHive.BaseRunnable):
    """subset a VCF by excluding/including (depending on the action value)
     the variants within the regions defined by a BED file"""

    def run(self):
        filepath = self.param('filepath')
        file = os.path.split(filepath)[1]
        work_dir = os.path.split(filepath)[0]

        if self.param_is_defined('work_dir'):
            work_dir = self.param('work_dir')

        self.warning('Analysing file: %s'% filepath)

        if self.param('verbose') == "True":
            print("Workdir is %s" % work_dir)

        vcfFilter = BCFTools(vcf=filepath, bcftools_folder=self.param('bcftools_folder'))
        vcffile = vcfFilter.subset_vcf(bed=self.param_required('bed'),
                                       outprefix=file+".exc.vcf.gz",
                                       outdir=work_dir,
                                       create_index=True,
                                       threads=self.param('threads'),
                                       action=self.param('action'))

        self.param('subset_file', vcffile)
        self.param('subset_file_ix', vcffile+'.tbi')

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow({'subset_file': self.param('subset_file')}, 1)
        self.dataflow({'subset_file_ix': self.param('subset_file_ix')}, 2)
