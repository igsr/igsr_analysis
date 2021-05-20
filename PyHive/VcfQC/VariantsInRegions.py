import eHive
import os
from VcfQC import VcfQC
from ReseqTrackDB import File
from ReseqTrackDB import ReseqTrackDB
from ReseqTrackDB import Attribute

class VariantsInRegions(eHive.BaseRunnable):
    """
    This runnable will get the number of variants in certain regions passed in a BED file

    Returns
    -------
    A param with the path to the file containing the coverage for the regions in the BED file
    """
    
    def run(self):

        filepath=self.param_required('filepath')

        self.warning('Analysing file: %s'% filepath)

        vcfQC = VcfQC(vcf=filepath,bedtools_folder=self.param_required('bedtools_folder'))

        outprefix=""

        if self.param_is_defined('work_dir'):
            file=os.path.split(filepath)[1]
            outprefix=self.param('work_dir')+"/"+file
        else:
            outprefix=filepath

        outfile=vcfQC.number_variants_in_region(region=self.param_required('region'),outprefix=outprefix)

        self.param('outfile', outfile)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'outfile' : self.param('outfile') }, 1)
