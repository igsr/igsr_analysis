import sys
import glob
import fnmatch
import eHive
from BamQC import BamQC
from ReseqTrackDB import *

class RunVerifyBamId(eHive.BaseRunnable):
    """run VerifyBAMID on a BAM file"""

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

        (path, filename) = os.path.split(filepath)

        population_id = filename.split(".")[3]

        if not population_id:
            raise Exception("Population could not be extracted from %s:" % filename)

        genotype_f = glob.glob(self.param_required('genotype_folder')+"/"+population_id+"*")

        if not len(genotype_f) == 1:
            self.warning("No population genotype file for %s" % filename)
            sys.exit()

        bam = BamQC(bam=filepath, verifybamid_folder=self.param_required('verifybamid_folder'))

        outfiles = bam.run_verifybamid(genotype_f[0], outprefix=filename,
                                       outdir=self.param_required('work_dir'))

        vfbamid_list = []

        for out in outfiles:
            index = []
            if fnmatch.fnmatch(out, '*.selfSM') or fnmatch.fnmatch(out, '*.selfRG'):
                with open(out) as f:
                    for line in f:
                        line = line.rstrip('\n')
                        if re.match(r"^#", line):
                            bits = line.split('\t')
                        else:
                            columns = line.split('\t')
                            sample = columns[bits.index('#SEQ_ID')]
                            rg = columns[bits.index('RG')]
                            freemix = columns[bits.index('FREEMIX')]
                            chipmix = columns[bits.index('CHIPMIX')]
                            #dropping trailing 0s
                            freemix = freemix.rstrip("0")
                            freemix = freemix.rstrip(".")
                            chipmix = (str(chipmix)[-2:] == '.0' and
                                       str(chipmix)[:-2] or str(chipmix))

                            name1 = "%s_%s:freemix" %(sample, rg)
                            name2 = "%s_%s:chipmix" %(sample, rg)
                            a1 = Attribute(table_name="file", other_id=fileO.dbID,
                                           name=name1, value=freemix)
                            a2 = Attribute(table_name="file", other_id=fileO.dbID,
                                           name=name2, value=chipmix)
                            vfbamid_list.append(a1.__dict__)
                            vfbamid_list.append(a2.__dict__)

        self.param('vfbamid_list', vfbamid_list)

    def write_output(self):
        self.warning('Work is done!')
        self.warning('{0} different Attributes were passed down'.
                     format(len(self.param('vfbamid_list'))))
        for attrb in self.param('vfbamid_list'):
            self.dataflow({'attrb': attrb}, 1)
