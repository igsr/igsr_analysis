import eHive

class SeedShapeit(eHive.BaseRunnable):
    """
    Class for seeding the PyHive::PipeConfig::INTEGRATION::Shapeit pipeline

    This pipeline requires a file with the following format:

    <input_bed.prefix>\t<chr>
    """
    def run(self):

        self.warning('Analysing file: %s'% self.param_required('filepath'))

        filepath = self.param_required('filepath')
        prefix = self.param_required('prefix')

        flist = [] # will store the list of tuples
        with open(filepath) as f:
            for line in f:
                line = line.rstrip('\n')
                flist.append((line.split('\t')[0], prefix+"."+line.split('\t')[1]))

        self.param('flist', flist)

    def write_output(self):
        self.warning('Work is done!')

        for f in self.param('flist'):
            self.dataflow({
                'input_bed' : f[0],
                'outprefix' : f[1]}, 1)
