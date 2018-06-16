import eHive
import subprocess
import os
import pdb
import sys
import random
import string

class VcfConcat(eHive.BaseRunnable):
    """Concat each of the VCF chunks into a single VCf"""
    
    def random_generator(self, size=6, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for x in range(size))
  
    def run(self):

        all_ixs = self.param_required('allixs')
        all_files = self.param_required('allchunks_files')

        if not os.path.isdir(self.param_required('work_dir')):
            os.makedirs(self.param_required('work_dir'))

        #log file
        if not os.path.isdir(self.param_required('log_dir')):
            os.makedirs(self.param_required('log_dir'))

        outprefix=os.path.split(self.param_required('outprefix'))[1]
        outprefix="{0}/{1}".format(self.param_required('work_dir'),outprefix)

        d = dict(zip(all_ixs, all_files))

        #write to log file
        logfile_url="{0}/merge_vcf.log".format(self.param_required('log_dir'))
        logf=open(logfile_url,'w')
        logf.write("#ix chunk\n")
        for k in sorted(d):
            logf.write("{0} {1}\n".format(k,d[k]))
        logf.close

        """Create tmp file for files to concat"""
        concat_file="%s/concat%s.txt"% (self.param_required('work_dir'),self.random_generator())

        f=open(concat_file,'w');
        for ix in sorted(d):
            f.write(d[ix]+"\n")
        f.close()    

        nt=1 #number of threads
        if self.param_is_defined('threads'):
            nt=self.param('threads')

        command="{0}/bcftools concat -f {1} --threads {2} ".format(self.param_required('bcftools_folder'),concat_file,nt)
        command= command+"-o {0} -O z".format(outprefix)

        if self.param('verbose')=="True":
            self.warning("Merging chunks")
            self.warning("Command used to merge chunks: %s" % command)
 
        try:
            subprocess.check_output(command,shell=True)
        except subprocess.CalledProcessError as e:
            self.warning("Something went wrong while merging the chunks")
            print(e.output)
            sys.exit(0)

        os.remove(concat_file)

        if self.param('verbose')=="True":
            self.warning("Merging the chunks: DONE")

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'merged_file' : self.param('outprefix') }, 1)




