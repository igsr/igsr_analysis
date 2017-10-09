import eHive
from BEDTools import BEDTools

class CoordFactory(eHive.BaseRunnable):
    """Run BEDTools 'make_windows' in order to create genomic chunks of a certain length (in bp)"""
    
    def run(self):
    
        bedtools_obj=BEDTools(bedtools_folder=self.param_required('bedtools_folder'))

        verbose=None
        if self.param_is_defined('verbose'):
            verbose=True
        else:
            verbose=False

        offset=None
        if self.param_is_defined('offset'):
            offset=self.param('offset')
            
        chrom=None
        if self.param_is_defined('chrom'):
            chrom=self.param('chrom')

        coord_list=bedtools_obj.make_windows(g=self.param_required('genome_file'), 
                                             w=self.param_required('window'), 
                                             s=offset,
                                             chrom=chrom,
                                             verbose=verbose)

        self.param('coord_list', coord_list)
       
    def write_output(self):
        self.warning('{0} chunks have been created'.format(len(self.param('coord_list'))))

        for chunk in self.param('coord_list'):
            self.dataflow( { 'chunk' : chunk }, 2)




