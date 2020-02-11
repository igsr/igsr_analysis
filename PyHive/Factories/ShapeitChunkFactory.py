import eHive
import glob
import os
import pdb

from VCF.VCFIntegration.Beagle import Beagle

class ShapeitChunkFactory(eHive.BaseRunnable):
    """Run makeBGLCHUNKS in order to create genomic chunks used to be run with Shapeit"""

    def run(self):

        vcf_i = Beagle(vcf=self.param_required('filepath'),
                       makeBGLCHUNKS_folder=self.param_required('makeBGLCHUNKS_folder'))

        outfile = ""
        if self.param_is_defined('work_dir'):
            if not os.path.isdir(self.param_required('work_dir')):
                os.makedirs(self.param_required('work_dir'))
            outfile += self.param('work_dir')+"/output.coords"

        #delete old Shapeit files
        for file in glob.glob("{0}/*shapeit*".format(self.param('work_dir'))):
            os.remove(file)

        verbose = None
        if self.param_is_defined('verbose'):
            verbose = True
        else:
            verbose = False

        window = 0
        if self.param_required('window') is not None:
            window = self.param('window')
        else:
            raise Exception("I need a 'window' arg")

        overlap = 0
        if self.param_required('overlap') is not None:
            overlap = self.param('overlap')
        else:
            raise Exception("I need an 'overlap' arg")

        vcf_i.make_beagle_chunks(window=self.param_required('window'),
                                 overlap=self.param_required('overlap'),
                                 outfile=outfile, verbose=verbose)

        chunks = []
        with open(outfile) as f:
            for line in f:
                line = line.rstrip('\n')
                # makeBGLCHUNKS does not work for chromosome X and will write
                # coordinates with an incorrect chromosome name ('0' in this case)
                if line.split('\t')[0] == '0':
                    coords = "{0}:{1}-{2}".format('X',
                                                  line.split('\t')[1],
                                                  line.split('\t')[2])
                else:
                    coords = "{0}:{1}-{2}".format(line.split('\t')[0],
                                                  line.split('\t')[1],
                                                  line.split('\t')[2])
                chunks.append(coords)

        self.param('chunks', chunks)

    def write_output(self):
        self.warning('{0} files have been created'.format(len(self.param('chunks'))))

        for chunk in self.param('chunks'):
            self.dataflow({'chunk': chunk}, 2)