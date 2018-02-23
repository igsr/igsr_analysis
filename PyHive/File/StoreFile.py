import eHive
import os
import pdb
import ast

from datetime import datetime
from ReseqTrackDB import File
from ReseqTrackDB import ReseqTrackDB


class StoreFile(eHive.BaseRunnable):
    """Store the file in the DB"""
    
    def param_defaults(self):
        return {
            'compression' : None
        }

    def run(self):
        
        self.warning('Storing file: %s'% self.param_required('filename'))

        # First, rename the file
        fileO=File(path=self.param_required('filename'),type=self.param_required('type'))
        filelayout=ast.literal_eval(self.param_required('filelayout'))
        newlayout=ast.literal_eval(self.param_required('newlayout'))

        compression=None
        if self.param('compression'):
            compression= self.param('compression')

        add_date=False
        if self.param('add_date'):
            add_date=True
    
        fileO.rename(filelayout= filelayout, newlayout= newlayout,
                     extension= self.param_required('extension'), add_date=add_date, 
                     compression= compression)

        new_f=File(path=self.param_required('filename'),type=self.param_required('type'))

        new_f.move(self.param('final_dir')+"/"+fileO.name)

        if self.param_is_defined('store'):
            hostname=self.param_required('hostname')
            username=self.param_required('username')
            db=self.param_required('db')
            port=self.param_required('port')
            pwd=self.param_required('pwd')

            reseqdb = ReseqTrackDB(host=hostname,user=username,port=port,pwd=pwd,db=db)

            new_f.store(reseqdb, do_md5= True, dry= False)

        self.param('stored_file', new_f.path)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'stored_file' : self.param('stored_file') }, 1)
