import eHive
import random
import os
import glob
import pdb
import string

def random_generator(size=6, chars=string.ascii_uppercase + string.digits, outprefix=None):
    '''
    Creates a random string having a size determined by the 'size' arg. If outprefix is not None
       then use a prefix with the random string

    Parameters
    ----------
    size : int, Required
           Size for the random string
    outprefix : str, Optional
                Add an outprefix for the returned random string
    '''
    astr=''.join(random.choice(chars) for x in range(size))
    if outprefix is not None:
        return "{0}_{1}".format(outprefix,astr)
    else:
        return astr

class ShortenFilePaths(eHive.BaseRunnable):
    '''
    Takes a list of files and shorten their paths by creating symbolic links

    Parameters
    ----------
    filelist : List, Required
               List with the path to the file/s that will be shortened
    work_dir : str, Required
               String to folder that will be used to put the files in

    Returns
    -------

    File with list of new shortened file paths
    '''

    def run(self):
        
        input_bams_str=""

        if not os.path.isdir(self.param_required('work_dir')):
            os.makedirs(self.param_required('work_dir'))
        
        filelist=self.param_required('filelist')

        short_files=[]

        for afile in filelist:

            folder_prefix="{0}short".format(random_generator(size=6))
            tmp_folder="{0}/{1}".format(self.param_required('work_dir'), folder_prefix)

            if not os.path.isdir(tmp_folder):
                os.makedirs(tmp_folder)

            # changing to tmp_folder dir
            os.chdir(tmp_folder);
        
            filenames="{0}/filenames.txt".format(tmp_folder)
            fw=open(filenames,'w')
            
            with open(afile) as f:
                for line in f:
                    old_id=line.rstrip("\n")
                    old_ix=old_id+".bai"
                    random_str=random_generator(size=6)
                    new_fileid="{0}.bam".format(random_str)
                    fw.write(new_fileid+"\n")
                    new_fileix="{0}.bam.bai".format(random_str)
                    os.symlink(old_id, new_fileid)
                    os.symlink(old_ix, new_fileix)
                
            fw.close()
            short_files.append(filenames)

        os.chdir(self.param_required('work_dir'))

        self.param('short_flist', short_files)

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow( { 'short_flist' : self.param('short_flist') }, 1)
