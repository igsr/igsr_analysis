import eHive
import os

from datetime import datetime
from ReseqTrackDB import File
from ReseqTrackDB import ReseqTrackDB


class StoreFile(eHive.BaseRunnable):
    """Store the file in the DB"""
    
    def param_defaults(self):
        return {
            'compression' : None
        }

    def fetch_input(self):
        hostname=self.param_required('hostname')
        username=self.param_required('username')
        db=self.param_required('db')
        port=self.param_required('port')
        pwd=self.param_required('pwd')

        reseqdb = ReseqTrackDB(host=hostname,user=username,port=port,pwd=pwd,db=db)

        self.param('reseqdb', reseqdb)

    def run(self):
        
        self.warning('Storing file: %s'% self.param_required('filename'))

        # First, rename the file
        newf=File(path=self.param_required('filename'),type=self.param_required('type'))
        oldlayot=self.param_required('oldlayout')
        newlayout=self.param_required('newlayout')

        compression=None
        if self.param('compression'):
            compression= self.param('compression')

        add_date=False
        if self.param('add_date'):
            add_date=True

        newf.rename(filelayout= self.param_required('oldlayout'), newlayout= self.param_required('newlayout'),
                    extension= self.param_required('extension'), add_date=add_date, compression= compression)

            
        self.param('stored_file', newf.path)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'stored_file' : self.param('stored_file') }, 1)
