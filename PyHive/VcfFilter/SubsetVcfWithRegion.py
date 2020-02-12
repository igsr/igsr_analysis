import eHive
import os
from VCFfilter.BCFTools import BCFTools

class SubsetVcfWithRegion(eHive.BaseRunnable):
    """subset a VCF by excluding/including (depending on the action value) the variants
    within the regions passed in as a parameter"""

    def run(self):
        filepath = self.param('filepath')
        file = os.path.split(filepath)[1]

        self.warning('Analysing file: %s'% filepath)

        work_dir = self.param_required('work_dir')
        vcfFilter = BCFTools(vcf=filepath, bcftools_folder=self.param('bcftools_folder'))

        region = self.param_required('region')
        self.warning('Subsetting file using region: %s'% region)
        if not os.path.isdir(work_dir):
            os.mkdir(work_dir)

        vcffile = vcfFilter.subset_vcf(region=region,
                                       outprefix=file,
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
