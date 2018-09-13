import eHive
import os
import pdb
import subprocess

class RunSamToolsMerge(eHive.BaseRunnable):
    """run Samtools merge on some BAM files"""
    
    def run(self):
        bamlist=self.param_required('bamlist')

        if len(bamlist)==1:
            self.param('merged_bam', bamlist[0])
        else:
            bamstr=" ".join(bamlist)
            region_str='_'.join(map(str,self.param_required('region')))
            
            merged_dir="{0}/merged_bams".format(self.param('work_dir'))
            if not os.path.isdir(merged_dir):
                os.makedirs(merged_dir)
        
            merged_f="{0}/merged_tpsed{1}.bam".format(merged_dir,region_str)
            cmd="{0}/samtools merge {1} {2}".format(self.param_required('samtools_folder'), merged_f, bamstr) 
            try:
                subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError as exc:
                print("Something went wrong while running Samtools Merge.\n"
                      "Command used was: %s" % cmd)
                raise Exception(exc.output)
            self.param('merged_bam', merged_f)
        
    def write_output(self):
        self.warning('Work is done!')

        self.dataflow( {'merged_bam' : self.param('merged_bam') }, 1)
