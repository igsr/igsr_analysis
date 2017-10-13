import eHive
import os
from VcfQC import VcfQC
from ReseqTrackDB import File
from ReseqTrackDB import ReseqTrackDB
from ReseqTrackDB import Attribute

class BcftoolsStats(eHive.BaseRunnable):
    """run BCFtools stats on a VCF file"""
    
    def param_defaults(self):
        return {
        }


    def __str_to_bool(self,s):
        if s == 'True':
            return True
        elif s == 'False':
            return False
        else:
            raise ValueError # evil ValueError!

    def fetch_input(self):
        hostname=self.param('hostname')
        username=self.param('username')
        db=self.param('db')
        port=self.param('port')
        pwd=self.param('pwd')

        reseqdb = ReseqTrackDB(host=hostname,user=username,port=port,pwd=pwd,db=db)

        self.param('reseqdb', reseqdb)

    def run(self):
        reseqdb=self.param('reseqdb')

        filepath=self.param_required('filepath')

        self.warning('Analysing file: %s'% filepath)

        file=os.path.split(filepath)[1]
        outfile=self.param_required('work_dir')+"/"+file

        vcfQC = VcfQC(vcf=filepath,bcftools_folder=self.param_required('bcftools_folder'))
    
        verbose=None
        if self.param_is_defined('verbose'):
            verbose=self.__str_to_bool(self.param('verbose'))

        params=dict()
        attr_suffix=""
        stats=None
        if self.param_is_defined('region'):
            attr_suffix="_{0}".format(self.param('region'))
            params['region']=self.param('region')
            
        if self.param_is_defined('filter_str'):
            attr_suffix += "_filt"
            params['filter_str']=self.param('filter_str')
        
        stats=vcfQC.stats(outpath=outfile,**params,verbose=verbose)

        if self.param_required('store_attributes')=='True':

            fileO=reseqdb.fetch_file_by_url(filepath)
            for attr,value in stats.summary_numbers.items():
                Attribute(table_name="file",other_id=fileO.dbID,name="STATS{0}_{1}".format(attr_suffix,attr),value=value).store(reseqdb)
            Attribute(table_name="file",other_id=fileO.dbID,name="STATS_ts_tv{0}".format(attr_suffix),value=stats.ts_tv).store(reseqdb)
            Attribute(table_name="file",other_id=fileO.dbID,name="STATS_ts_tv_1stalt{0}".format(attr_suffix),value=stats.ts_tv_1stalt).store(reseqdb)
            Attribute(table_name="file",other_id=fileO.dbID,name="STATS_no_singleton_snps{0}".format(attr_suffix),value=stats.no_singleton_snps).store(reseqdb)

        self.param('stats_file', stats.filename)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'stats_file' : self.param('stats_file') }, 1)
