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

        #store file
        newf=File(path=self.param_required('filename'),type=self.param_required('type'))
        layout_dict=self.param_required('layout_dict')
        newlayout=self.param_required('newlayout')

        newpath=self.param_required('final_dir')+"/"
        bits=[ layout_dict[i] for i in newlayout]

        if self.param_required('add_date')=='True':
            now = datetime.now()
            bits.append(now.strftime('%Y%m%d') )
        
        bits.append(self.param_required('extension'))

        if self.param('compression'):
            bits.append(self.param('compression'))

        newpath+=".".join(bits)        
        newf.move(self.param('reseqdb'),do_md5=True,newpath=newpath)
            
        self.param('stored_file', newf.path)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'stored_file' : self.param('stored_file') }, 1)
