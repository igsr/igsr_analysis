import eHive
import random
import os
import glob
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
    astr = ''.join(random.choice(chars) for x in range(size))
    if outprefix is not None:
        return "{0}_{1}".format(outprefix, astr)
    else:
        return astr

class ShortenFiles(eHive.BaseRunnable):
    '''
    Takes a list of files and shorten their paths by creating symbolic links

    Parameters
    ----------
    filelist : str, Required
               Path to files that will be shortened
    work_dir : str, Required
               String to folder that will be used to put the files in

    Returns
    -------

    File with list of new shortened files
    '''

    def run(self):

        tmp_folder = self.param_required('work_dir')+"/tmp_sh"
        if not os.path.isdir(self.param_required('work_dir')):
            os.makedirs(self.param_required('work_dir'))
        if not os.path.isdir(tmp_folder):
            os.makedirs(tmp_folder)

        #delete files from previous runs
        for file in glob.glob(tmp_folder+"/*"):
            os.unlink(file)

        # changing to tmp_folder dir
        os.chdir(tmp_folder)

        fw = open(tmp_folder+"/files_to_shorten.txt", 'w')

        with open(self.param_required('filelist')) as f:
            for line in f:
                old_id = line.rstrip("\n")
                old_ix = old_id+".bai"
                random_str = random_generator(size=6)
                new_fileid = "{0}.bam".format(random_str)
                fw.write(new_fileid+"\n")
                new_fileix = "{0}.bam.bai".format(random_str)
                os.symlink(old_id, new_fileid)
                os.symlink(old_ix, new_fileix)

        fw.close()

        os.chdir(self.param_required('work_dir'))

        self.param('short_f', tmp_folder+"/files_to_shorten.txt")

    def write_output(self):
        self.warning('Work is done!')
        self.dataflow({'short_f': self.param('short_f')}, 1)
