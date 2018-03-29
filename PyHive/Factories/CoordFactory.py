import eHive
import pdb
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

        subtract=None
        if self.param_is_defined('subtract'):
            subtract=self.param('subtract')

        rextend=None
        if self.param_is_defined('rextend'):
            rextend=self.param('rextend')
          
        coord_list=bedtools_obj.make_windows(g=self.param_required('genome_file'), 
                                             w=self.param_required('window'), 
                                             s=offset,
                                             subtract=subtract,
                                             rextend=rextend,
                                             verbose=verbose)
        chunks=[]

        ix=1
        for c in coord_list:
            if self.param_is_defined('ix'):
                if ix==self.param('ix'):
                    chunks.append(
                        {
                            'chunk': c,
                            'ix': ix
                        })
            else:
                chunks.append(
                    {
                        'chunk': c,
                        'ix': ix
                    })
            ix+=1

        self.param('chunks', chunks)
       
    def write_output(self):
        self.warning('{0} chunks have been created'.format(len(self.param('chunks'))))

        self.dataflow(self.param('chunks'), 2)





