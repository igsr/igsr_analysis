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
        
        filenames=self.param_required('short_flist')

        bamlist=[]

        for afile in filenames:
        
            afile_dir=os.path.dirname(os.path.abspath(afile))

            bam_dir="{0}/tposed_bams".format(afile_dir)

            if not os.path.isdir(bam_dir):
                os.makedirs(bam_dir)

            #change to folder containing shortened BAMs
            os.chdir(afile_dir)

            region="{0}:{1}-{2}".format(self.param_required('region')[0], self.param_required('region')[1], self.param_required('region')[2])
            region_str='_'.join(map(str,self.param_required('region')))
            output_file="{0}/{1}.{2}.tranposed.bam".format(bam_dir, self.param_required('outprefix'), region_str)

            cmd="cat {0} | xargs -s 2000000 {1}/transpose_bam -r {2} -i -o {3}".format(
                afile,
                self.param_required('transposebam_folder'),
                region,
                output_file)

            try:
                subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError as e:
                print(e.output)
                sys.exit(0)

            bamlist.append(output_file)

        #change back to working dir
        os.chdir(self.param_required('work_dir'))

        self.param('out_bamlist', bamlist)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'out_bamlist' : self.param('out_bamlist') }, 1)

        
