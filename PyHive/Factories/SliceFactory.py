import eHive
import os
from Coord import *

class SliceFactory(eHive.BaseRunnable):
    """Slice a file into windows of a certain length"""

    def param_defaults(self):
        return {
        }

    def run(self):

        self.warning('Creating slices for sequences in: %s'% self.param_required('faix'))

        if os.path.isfile(self.param_required('faix')) is False:
            raise Exception("File does not exist")

        slices = []
        with open(self.param_required('faix')) as f:
            for line in f:
                elms = line.split("\t")
                c = Coord(id=elms[0], start=0, end=int(elms[1]))
                chrlist = c.make_windows(step=self.param_required('slice_size'))
                for i in chrlist:
                    slices.append({'id':i.id, 'start':i.start, 'end':i.end})

        self.param('slices', slices)

    def write_output(self):
        self.warning('{0} chros jobs have been created'.format(len(self.param('slices'))))
        self.dataflow(self.param('slices'), 2)
