import eHive
import os

class SeedFile(eHive.BaseRunnable):
    """Class for seeding a pipeline with the files listed in 'filepath'"""

    def run(self):

        self.warning('Analysing file: %s'% self.param_required('filepath'))
        
        filepath=self.param_required('filepath')
        
        flist=[] # will store the list of files
        with open(filepath) as f:
            for line in f:
                line=line.rstrip('\n')
                flist.append(line)
        
        self.param('flist', flist)

    def write_output(self):
        self.warning('Work is done!')

        for f in self.param('flist'):
            self.dataflow( { 'filepath' : f }, 1)
