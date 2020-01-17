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

        if len(all_ixs)!= len(all_files):
            raise Exception("Number of indexes does not match with number files")

        data={}
        for i in all_files:
            coord_str=os.path.basename(i).split('.')[-3]
            chr,start,end=coord_str.split('_')
            if not chr in data:
                data[chr]=[(int(start),int(end),i)]
            else:
                data[chr].append((int(start),int(end),i))

        sorted_files=[]
        for i in ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
                  'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
                  'chr18','chr19','chr20','chr21','chr22','chrX']:
            if i in data:
                coords=data[i]
                s_list=sorted(coords, key=lambda element: (element[1], element[2]))
                for f in s_list:
                    sorted_files.append(f[2])

        if not os.path.isdir(self.param_required('work_dir')):
            os.makedirs(self.param_required('work_dir'))

        outprefix=os.path.split(self.param_required('outprefix'))[1]
        outprefix="{0}/{1}".format(self.param_required('work_dir'),outprefix)

        """Create tmp file for files to concat"""
        concat_file="%s/concat%s.txt"% (self.param_required('work_dir'),self.random_generator())
        f=open(concat_file,'w');
        for k in sorted_files:
            f.write(k+"\n")
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




