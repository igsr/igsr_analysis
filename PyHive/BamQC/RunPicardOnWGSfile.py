import eHive
from BamQC import BamQC
from ReseqTrackDB import *

class RunPicardOnWGSfile(eHive.BaseRunnable):
    """run Picard's CollectWgsMetrics on a WGS BAM file"""

    def param_defaults(self):
        return {
        }

    def fetch_input(self):
        hostname = self.param('hostname')
        username = self.param('username')
        db = self.param('db')
        port = self.param('port')
        pwd = self.param('pwd')

        reseqdb = ReseqTrackDB(host=hostname,user=username,port=port,pwd=pwd,db=db)

        self.param('reseqdb', reseqdb)

    def run(self):
        filepath=self.param_required('filepath')
        reseqdb=self.param('reseqdb')

        fileO=reseqdb.fetch_file_by_url(filepath)

        self.warning('Analysing file: %s'% filepath)
        
        (path,filename)=os.path.split(filepath)
        
        bam = BamQC(bam=filepath,java_folder=self.param_required('java_folder'),picard_folder=self.param_required('picard_folder'))

        CHsMetricsO=bam.run_CollectWgsMetrics(self.param_required('reference'),outfile=self.param_required('work_dir')+"/"+filename+".CollectWgsMetrics.metrics.txt")
        #set the X-axis limit depending on the coverage
        xlims=None
        if (float(CHsMetricsO.metrics.get('MEDIAN_COVERAGE'))<20):
            xlims=[0,50]
        elif (float(CHsMetricsO.metrics.get('MEDIAN_COVERAGE'))>=20):
            xlims=[0,100]
        
        picard_list=[]
        for k,v in CHsMetricsO.metrics.items(): 
            attrb=Attribute(table_name="file",other_id=fileO.dbID,name="PICARD:"+k,value=v)
            picard_list.append(attrb.__dict__)

        CHsMetricsO.create_cov_barplot(self.param_required('work_dir')+"/"+filename+".CollectWgsMetrics.barplot.pdf",xlim=xlims)
                            
        self.param('picard_list', picard_list)
        
    def write_output(self):
        self.warning('Work is done!')
        self.warning('{0} different Attributes were passed down'.format(len(self.param('picard_list'))))
        for attrb in self.param('picard_list'):
            self.dataflow({ 'attrb' : attrb}, 1)
