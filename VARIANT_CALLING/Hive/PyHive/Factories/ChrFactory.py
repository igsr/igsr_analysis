import eHive


class ChrFactory(eHive.BaseRunnable):
    '''
    Factory that creates 1 job per chromosome found in a Fasta index file.
    
    Returns
    -------
    This factory returns a list of dictionaries, one dict per chromosome
    '''

    def fetch_input(self):
        faix = self.param_required('faix')

        files=[]
        ix=1
        for line in open(faix):
            if line.startswith("\n"):
                continue
            chro=line.split('\t')[0]
            files.append(
                {
                    'chr': chro,
                    'ix': ix
                }
            )
            ix+=1
        self.param('files', files)
        
    def write_output(self):
        self.warning('{0} chros jobs have been created'.format(len(self.param('files'))))
        self.dataflow(self.param('files'), 2)

