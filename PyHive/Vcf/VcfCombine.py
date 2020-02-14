import eHive

from VCF.VcfUtils import VcfUtils

class VcfCombine(eHive.BaseRunnable):
    '''
    Combine different VCFs generated by different callers into a single VCF

    This runnable will call run GATK CombineVariants with the following options by default:

    -env: Exclude sites where no variation is present after merging
    -sites_only: Drop genotype information and consider only the sites
    --filteredAreUncalled: this flag causes filtered variants (i.e. variant
    records where the FILTER field is populated by
    something other than PASS or a dot) to be omitted from the output.
    --genotypemergeoption='UNIQUIFY'  Make all sample genotypes unique by file.
    '''

    def run(self):

        vcf_paths = self.param_required('allfiles2combine')
        labels = self.param_required('alldatasets2combine')

        vcf_utils = VcfUtils(vcflist=vcf_paths,
                             bcftools_folder=self.param_required('bcftools_folder'),
                             gatk_folder=self.param_required('gatk_folder'),
                             tmp_dir=self.param('tmp_dir'))

        ginterval = None
        if self.param_is_defined('ginterval'):
            ginterval = self.param('ginterval')

        outfile = vcf_utils.combine(labels=labels, reference=self.param_required('reference'),
                                    threads=self.param_required('threads'),
                                    outprefix=self.param_required('outprefix'),
                                    outdir=self.param_required('work_dir'),
                                    compress=True,
                                    genotypemergeoption='UNIQUIFY',
                                    ginterval=ginterval,
                                    options=['-env', '-sites_only', '--filteredAreUncalled'])
        self.param('out_vcf', outfile)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow({'out_vcf': self.param('out_vcf')}, 1)
