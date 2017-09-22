import eHive
from VCFIntegration.Beagle import Beagle

class BeagleChunkFactory(eHive.BaseRunnable):
    """Run makeBGLCHUNKS in order to create genomic chunks used to be run with Beagle"""
    
    def run(self):
    
        vcf_i=Beagle(vcf=self.param_required('filepath'), makeBGLCHUNKS_folder=self.param_required('makeBGLCHUNKS_folder'))

        outfile=""

        if self.param_is_defined('work_dir'):
            outfile+=self.param('work_dir')+"/output.coords"

        verbose=None
        if self.param_is_defined('verbose'):
            verbose=True
        else:
            verbose=False

        if self.param_is_defined('correct'):
            outfile=vcf_i.make_beagle_chunks(window=self.param_required('window'),overlap=self.param_required('overlap'),outfile=outfile,correct=True,
                                     chrname=self.param('chro'),verbose=verbose)
        else:
            vcf_i.make_beagle_chunks(window=self.param_required('window'),overlap=self.param_required('overlap'),outfile=outfile,correct=False,
                                     verbose=verbose)

        chunks=[]
        with open(outfile) as f:
            for line in f:
                line=line.rstrip('\n')
                coords="{0}:{1}-{2}".format(line.split('\t')[0],
                                            line.split('\t')[1],
                                            line.split('\t')[2])
                chunks.append(coords)

        self.param('chunks', chunks)
       
    def write_output(self):
        self.warning('{0} files have been created'.format(len(self.param('chunks'))))

        for chunk in self.param('chunks'):
            self.dataflow( { 'chunk' : chunk }, 1)




