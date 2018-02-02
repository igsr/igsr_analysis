import eHive
import os
import pdb
import ast

from ReseqTrackDB import File
from ReseqTrackDB import ReseqTrackDB


class SplitFile(eHive.BaseRunnable):
    """Split File into the different bits and map an annotation to each bit"""
    
    def param_defaults(self):
        return {
        }

    def fetch_input(self):
        filename = self.param_required('filename')

        hostname=self.param_required('hostname')
        username=self.param_required('username')
        db=self.param_required('db')
        port=self.param_required('port')
        pwd=self.param_required('pwd')

        reseqdb = ReseqTrackDB(host=hostname,user=username,port=port,pwd=pwd,db=db)

        fileO=reseqdb.fetch_file_by_filename(filename)
        self.param('file_object', fileO)

    def run(self):
        
        filename=self.param('file_object').name
        filepath=self.param('file_object').path
        self.warning('Splitting file: %s'% filename)

        bits=filename.split('.')
        file_layout=ast.literal_eval(self.param_required('filelayout'))
    
        if len(bits)!=len(file_layout):
            print("Passed file contains the following bits: {0}".format(",".join(bits)))
            print("Specified layout contains the following bits: {0}".format(",".join(file_layout)))
            raise Exception("Length is not the same for filename bits and its associated annotations that are passed using the file_layout param")
        
        d = dict(zip(file_layout, bits))
        self.param('layout_dict',d)
        self.param('filepath', filepath)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'layout_dict' : self.param('layout_dict') }, 1)
        self.dataflow( { 'filepath' : self.param('filepath') }, 1)
