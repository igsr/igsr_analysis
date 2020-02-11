import eHive

class SeedVCFIntegration(eHive.BaseRunnable):
    """Class for seeding the VCFIntegration pipeline"""

    def run(self):

        self.warning('Analysing file: %s'% self.param_required('filepath'))

        filepath = self.param_required('filepath')

        flist = [] # will store the list of tuples
        with open(filepath) as f:
            for line in f:
                line = line.rstrip('\n')
                if len(line.split('\t')) != 2:
                    print(line.split('\t'))
                    print("Error getting the right label/file pair for line:\n{0}".format(line))
                    print("Check if there is any extra tab space between label/file\n")
                    raise Exception()
                flist.append({
                    'dataset': line.split('\t')[0],
                    'file': line.split('\t')[1]
                })

        self.param('flist', flist)

    def write_output(self):

        self.warning('Work is done!')
        self.dataflow(self.param('flist'), 2)
