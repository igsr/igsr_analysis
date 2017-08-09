import eHive
import os
import subprocess

from datetime import datetime
from ReseqTrackDB import File
from ReseqTrackDB import ReseqTrackDB


class CompressBybgzip(eHive.BaseRunnable):
    """Store the file in the DB"""
    
    def param_defaults(self):
        return {
            'compression' : None
        }

    def run(self):
        
        self.warning('Compressing file: %s'% self.param_required('filename'))

        filename=self.param_required('filename')

        basename=os.path.basename(filename).split('.')[0]
        outfile=filename+".bgz"

        command = "{0}/bgzip -c {1} > {2}".format(self.param_required('bgzip_folder'), filename, outfile)

        try:
            subprocess.check_output(command,shell=True)
        except subprocess.CalledProcessError as e:
            print(e.output)

        if self.param_required('delete_uncompressed')=='True':
            os.remove(filename)

        if self.param_required('create_index')=='True':
            outfile_ix=outfile+".tbi"
            command_ix="{0}/tabix {1}".format(self.param_required('tabix_folder'),outfile)
            try:
                subprocess.check_output(command_ix,shell=True)
                self.param('compressed_file_ix', outfile_ix)
            except subprocess.CalledProcessError as e:
                print(e.output)
            
        
        self.param('compressed_file', outfile)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'compressed_file' : self.param('compressed_file') }, 1)
        self.dataflow( { 'compressed_file_ix' : self.param('compressed_file_ix') }, 1)
