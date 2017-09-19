import eHive
import tempfile

class CollectRawF(eHive.BaseRunnable):
    """Print into a single file the different *.raw files corresponding to each of the samples """

    def run(self):
        all_files = self.param_required('allraws_files')
        
        outfile=tempfile.NamedTemporaryFile(dir=self.param_required('work_dir'),delete=False,prefix='combined_rawlist')

        f=open(outfile.name,'w');
        for raw_f in all_files:
            f.write(raw_f+"\n")

        #add a newline at the end of the file
        f.write("\n")
        f.close
        
        self.param('outfile', outfile.name)


    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'outfile' : self.param('outfile') }, 1)
