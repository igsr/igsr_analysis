import eHive
import os
import pdb
import subprocess

class RunSamToolsIndex(eHive.BaseRunnable):
    """run SamTools index on a BAM file"""
    
    def run(self):
        bamfile=self.param_required('bamfile')
        
        cmd="{0}/samtools index {1}".format(self.param_required('samtools_folder'), bamfile) 
        
        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as exc:
            print("Something went wrong while running Samtools Merge.\n"
                  "Command used was: %s" % cmd)
            raise Exception(exc.output)
        
        self.param('bam_ix', bamfile+".bai")
            
    def write_output(self):
        self.warning('Work is done!')

        self.dataflow( {'bam_ix' : self.param('bam_ix') }, 1)
