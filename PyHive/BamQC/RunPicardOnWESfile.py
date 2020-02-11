import eHive
from BamQC import BamQC
from ReseqTrackDB import *

class RunPicardOnWESfile(eHive.BaseRunnable):
    """run Picard's CollectHsMetrics on a WES BAM file"""

    def fetch_input(self):
        hostname = self.param('hostname')
        username = self.param('username')
        db = self.param('db')
        port = self.param('port')
        pwd = self.param('pwd')

        reseqdb = ReseqTrackDB(host=hostname, user=username, port=port, pwd=pwd, db=db)

        self.param('reseqdb', reseqdb)

    def run(self):
        filepath = self.param_required('filepath')
        reseqdb = self.param('reseqdb')

        fileO = reseqdb.fetch_file_by_url(filepath)

        self.warning('Analysing file: %s'% filepath)

        filename = os.path.split(filepath)[1]

        bam = BamQC(bam=filepath, java_folder=self.param_required('java_folder'),
                    picard_folder=self.param_required('picard_folder'))

        CHsMetricsO = bam.run_CollectHsMetrics(baits_file=self.param_required('baits_file'),
                                               outfile=self.param_required('work_dir')
                                               +"/"+filename+".CollectHsMetrics.metrics.txt")

        picard_list = []
        for k, v in CHsMetricsO.metrics.items():
            if v == "":
                continue
            attrb = Attribute(table_name="file", other_id=fileO.dbID, name="PICARD:"+k, value=v)
            picard_list.append(attrb.__dict__)

        CHsMetricsO.create_cov_barplot(self.param_required('work_dir')+
                                       "/"+filename+".CollectHsMetrics.barplot.pdf")

        self.param('picard_list', picard_list)

    def write_output(self):
        self.warning('Work is done!')
        self.warning('{0} different Attributes were passed down'.
                     format(len(self.param('picard_list'))))
        for attrb in self.param('picard_list'):
            self.dataflow({'attrb': attrb}, 1)
