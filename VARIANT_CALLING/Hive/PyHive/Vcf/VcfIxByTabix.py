import eHive
import subprocess
import os
import sys

class VcfIxByTabix(eHive.BaseRunnable):
    """Create an index for a VCF file using Tabix"""
  
    def run(self):

        outfile=self.param_required('filepath')+".tbi"
        command="{0}/tabix {1}".format(self.param_required('tabix_folder'),self.param_required('filepath'))
 
        try:
            subprocess.check_output(command,shell=True)
        except subprocess.CalledProcessError as e:
            self.warning("Something went wrong while creating the index")
            print(e.output)
            sys.exit(0)

        self.param('outfile', outfile)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'vcf_ix' : self.param('outfile') }, 1)




