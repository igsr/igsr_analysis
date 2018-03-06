import eHive
import os
import pdb
import re
import subprocess

class TransposeBam(eHive.BaseRunnable):
    '''
    Takes BAM files from a list and creates transposed BAMs for a certain genomic region

    Parameters
    ----------
    filelist : str, Required
               Path to list of shortened bam files generated with ShortenFiles
    transposebam_folder : str, Required
                          Path to folder containing the transpose_bam binary
    region : str, Required
             Genomic region that will be transposed. i.e. chr1:1000-2000
    outprefix : str, Required
                Prefix used for output transposed bam file
    work_dir : str, Required
               String to folder that will be used to put the files in
    '''

    def run(self):
        
        input_bams_str=""

        if not os.path.isdir(self.param_required('work_dir')):
            os.makedirs(self.param_required('work_dir'))

        region="{0}:{1}-{2}".format(self.param_required('region')[0], self.param_required('region')[1], self.param_required('region')[2])
        region_str='_'.join(map(str,self.param_required('region')))
        output_file="{0}/{1}.{2}.tranposed.bam".format(self.param_required('work_dir'), self.param_required('outprefix'), region_str)

        cmd="cat {0} | xargs -s 2000000 {1}/transpose_bam -r {2} -i -o {3}".format(
            self.param_required('filelist'),
            self.param_required('transposebam_folder'),
            region,
            output_file)

        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            print(e.output)
            sys.exit(0)

        self.param('out_bam', output_file)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'out_bam' : self.param('out_bam') }, 1)

        
