import eHive
import os

from VariantCalling import BCFTools

class BCFTools_caller(eHive.BaseRunnable):
    '''
    Run BCFTools mpileup/call on a BAM file

    eHive runnable for BCFTools mpileup/call

    Parameters
    ----------
    work_dir : str, Required
         Path to the working directory
    outprefix : str, Required
         String used as prefix for the output file
    annots : str, Required
             mpileup parameter
             String representing a comma separated list of
             annotations used to decorate the VCF
    chunk: str, Required
         Interval used by BCFTools for the variant calling
    bamlist: str, Required
         Path to BAM file used for the variant calling
    reference: str, Required
         Path to the Fasta file used to generate the BAM alignment file
    bcftools_folder: str, Required
                     Path to folder containing the BCFTools binary
    E: Boolean, Optional
       mpileup parameter
       Recalculate BAQ on the fly, ignore existing BQ tags. Default=False
    F : float, Optional
        mpileup parameter
        Minimum fraction of gapped reads. Default=0.002
    p: Boolean, Optional
       mpileup parameter
       Apply -m and -F thresholds per sample to increase sensitivity of calling.
       By default both options are applied to reads pooled from all samples
    d : int, Optional
        mpileup parameter
        At a position, read maximally INT reads per input file. Default=250
    C : int, Optional
        mpileup parameter
        Coefficient for downgrading mapping quality for reads containing
        excessive mismatches. Default=50
    P : str, Optional
        mpileup parameter
        Comma-delimited list of patforms (determined by @RG-PL) from which indel
        candidates are obtained. Default= ILLUMINA
    m_pileup: int, Optional
              mpileup parameter
              Minimum number gapped reads for indel candidates. Default=1
    m_call : boolean, Optional
             call parameter
             alternative modelfor multiallelic and rare-variant calling designed
             to overcome known limitations in -c calling model
    v : bool, Optional
        call parameter
        output variant sites only
        Default: False
    S : str, Optional
        call parameter
        File of sample names to include or exclude if prefixed with "^"
    threads: int, Optional
         Number of CPUs used by the caller
         Default=1
    verbose : str, Optional
              Print command line. Possible values are 'True' or 'False'

    Returns
    -------

    Path to VCF file
    '''

    def run(self):

        if not os.path.isdir(self.param_required('work_dir')):
            os.makedirs(self.param_required('work_dir'))

        outprefix = os.path.split(self.param_required('outprefix'))[1]

        chrom = self.param_required('chunk')[0]
        start = self.param_required('chunk')[1]
        end = self.param_required('chunk')[2]

        region = "{0}:{1}-{2}".format(chrom, start, end)

        outfile = "{0}/{1}.bcftools".format(self.param_required('work_dir'),
                                            outprefix)

        bcftools_object = BCFTools(bam=self.param_required('bam'),
                                   reference=self.param_required('reference'),
                                   bcftools_folder=self.param_required('bcftools_folder'))

        nt = 1
        if self.param_is_defined('threads'):
            nt = self.param('threads')

        E = False
        if self.param_is_defined('E'):
            E = True

        p = False
        if self.param_is_defined('p'):
            p = True

        m_call = False
        if self.param_is_defined('m_call'):
            m_call = True

        m_pileup = 1
        if self.param_is_defined('m_pileup'):
            m_pileup = self.param('m_pileup')

        F = 0.002
        if self.param_is_defined('F'):
            F = self.param('F')

        d = 250
        if self.param_is_defined('d'):
            d = self.param('d')

        C = 50
        if self.param_is_defined('C'):
            C = self.param('C')

        P = "ILLUMINA"
        if self.param_is_defined('P'):
            P = self.param('P')

        S = None
        if self.param_is_defined('S'):
            S = self.param('S')

        v = False
        if self.param_is_defined('v'):
            v = True

        verbose = None
        if self.param_is_defined('verbose'):
            verbose = True
        else:
            verbose = False

        outfile = bcftools_object.run_bcftools(outprefix=outfile,
                                               annots=self.param_required('annots'),
                                               E=E,
                                               p=p,
                                               F=F,
                                               d=d,
                                               C=C,
                                               P=P,
                                               S=S,
                                               m_pileup=m_pileup,
                                               m_call=m_call,
                                               r=region,
                                               v=v,
                                               verbose=verbose)

        self.param('out_vcf', outfile)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow({'out_vcf': self.param('out_vcf')}, 1)
