import eHive
import os
import pdb
import ast

from ReseqTrackDB import File
from ReseqTrackDB import ReseqTrackDB


class SplitFile(eHive.BaseRunnable):
    """Split File into the different bits and map an annotation to each bit"""
    
    def fetch_input(self):
        filename = self.param_required('filename')

        fileO=None

        if self.param_is_defined('hostname') and self.param_is_defined('username'):
            hostname=self.param('hostname')
            username=self.param('username')
            db=self.param('db')
            port=self.param('port')
            pwd=self.param('pwd')

            reseqdb = ReseqTrackDB(host=hostname,
                                   user=username,
                                   port=port,
                                   pwd=pwd,
                                   db=db)

            fileO=reseqdb.fetch_file_by_filename(filename)
        else:
            fileO=File(
                path=filename, 
                type='PHASED_VCF',
            )
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
        self.dataflow( { 
            'layout_dict' : self.param('layout_dict'),
            'filepath' : self.param('filepath')
        }, 1)

