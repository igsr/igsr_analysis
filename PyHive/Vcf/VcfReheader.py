import eHive
import os
import pdb
import tempfile
from VcfUtils import VcfUtils

class VcfReheader(eHive.BaseRunnable):
    """Modify the header in a VCF file"""

    def run(self):
        pdb.set_trace()
        filepath=self.param_required('filepath')
        newheader=self.param_required('newheader')

        self.warning('Analysing file: %s'% filepath)

        basename=os.path.split(filepath)[1]
        work_dir=None
        if self.param_is_defined('work_dir'):
            work_dir=self.param_required('work_dir')
        else:
            work_dir=os.path.split(filepath)[0]

        vcf_object=VcfUtils(vcf=filepath, bcftools_folder=self.param_required('bcftools_folder'))

        outprefix="{0}/{1}".format(work_dir, basename)

        outfile=""

        if self.param_is_defined('samplename'):
            samplename=self.param('samplename')

            try:
                tfile = tempfile.NamedTemporaryFile(mode='w',dir='data/')
                tfile.write(samplename+"\n")
                tfile.flush()
                outfile=vcf_object.reheader(newheader=self.param_required('newheader'), samplefile=tfile.name,
                                            outprefix=outprefix)
            finally:
                # Automatically cleans up the file
                tfile.close()
        elif self.param_is_defined('samplefile'):
            outfile=vcf_object.reheader(newheader=self.param_required('newheader'), samplefile=self.param('samplefile'),
                                        outprefix=outprefix)
        else:
            try:
                outfile=vcf_object.reheader(newheader=self.param_required('newheader'), outprefix=outprefix)
            except Exception as error:
                print('Caught this error: ' + repr(error))
                

        self.param('vcf_f', outfile)

    def write_output(self):
        self.warning('Work is done!')

        self.dataflow( {'vcf_f' : self.param('vcf_f') }, 1)
