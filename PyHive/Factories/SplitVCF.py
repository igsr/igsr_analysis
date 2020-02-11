import eHive
import os
import sys
import subprocess
import string
import random
import glob
from Coord import *

def random_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for x in range(size))

class SplitVCF(eHive.BaseRunnable):
    """Split a VCF into chunks containing a certain number of lines"""

    def run(self):

        self.warning('Split file: %s'% self.param_required('filepath'))

        ifile = self.param_required('filepath')
        outdir = self.param_required('work_dir')
        number_of_lines = self.param_required('number_of_lines')
        bgzip_folder = self.param_required('bgzip_folder')

        ## getting the header

        cmd1 = "zcat {0} | head -n 10000 |grep '^#' | " \
               "bgzip > {1}/header.txt.gz".format(ifile, outdir)

        if self.param('verbose') == "True":
            self.warning("Getting the header")
            self.warning("Command used to get the header: %s" % cmd1)

        try:
            subprocess.check_output(cmd1, shell=True)
        except subprocess.CalledProcessError as e:
            self.warning("Something went wrong while trying to the get the header")
            print(e.output)
            sys.exit(0)

        if self.param('verbose') == "True":
            self.warning("Getting the header: DONE")

        outprefix = outdir+"/"+random_generator(size=8)

        ## splitting

        cmd2 = "zgrep -v '^#' {0} | split -l {2} - --filter='{3}/bgzip > {1}$FILE.tmp.vcf.gz'".\
            format(ifile, outprefix, number_of_lines, bgzip_folder)

        if self.param('verbose') == "True":
            self.warning("Splitting the file")
            self.warning("Command used to split the file: %s" % cmd2)

        try:
            subprocess.check_output(cmd2, shell=True)
        except subprocess.CalledProcessError as e:
            self.warning("Something went wrong while splitting the file")
            print(e.output)
            sys.exit(0)

        if self.param('verbose') == "True":
            self.warning("Splitting the file: DONE")

        files = []

        #adding the header to each chunk
        for file in sorted(glob.glob(outprefix+"*")):
            final_file = file.replace('.tmp', '')
            cmd3 = "zcat {0}/header.txt.gz {1} | {3}/bgzip -c > {2}".\
                format(outdir, file, final_file, bgzip_folder)
            if self.param('verbose') == "True":
                self.warning("Adding the header: Command used is %s" % cmd3)
            try:
                subprocess.check_output(cmd3, shell=True)
            except subprocess.CalledProcessError as e:
                self.warning("Something went wrong while adding the header to each chunk file")
                print(e.output)
                sys.exit(0)
            os.remove(file)
            files.append({'filepath': final_file})

        self.param('files', files)

    def write_output(self):
        self.warning('{0} chunks have been created'.format(len(self.param('files'))))

        if self.param('verbose') == "True":
            for f in self.param('files'):
                self.warning("Split file is %s" % f)
        self.dataflow(self.param('files'), 2)
