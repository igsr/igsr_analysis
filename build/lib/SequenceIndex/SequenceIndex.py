'''
Created on 28 Nov 2016

@author: ernesto
'''

import re
from collections import defaultdict

class SequenceIndex(object):
    '''
    Class representing an index file as it is represented in the IGSR project
    '''

    def __init__(self, filepath):
        '''
        Constructor

         Class variables
        ---------------
        filepath : str, required
                Path to index file
        columns : list
                List with column names in the file
        '''

        self.filepath = filepath

        with open(self.filepath) as f:
            for line in f:
                line = line.rstrip('\n')
                if re.match(r"^#{1}\w+", line):
                    bits = line.split('\t')
                    self.columns = bits
                    break

    def runs_per_sample(self, analysis_group=None):
        '''
        This function will return a dictionary for which each key will be a different
        sample and the values will be a set with run ids for that sample.

        IMPORTANT: Only index lines corresponding to the ILLUMINA instrument
        platform will be considered

        Parameters
        ----------
        analysis_group : str, optional
                         Only retrieve the run_ids for this particular
                         analysis_group (i.e. 'low coverage')

        '''

        data = defaultdict(set)

        with open(self.filepath) as f:
            for line in f:
                line = line.rstrip('\n')
                if not line.startswith('#'):
                    if not line.split('\t')[self.columns.index('INSTRUMENT_PLATFORM')] \
                    == 'ILLUMINA': continue
                    if analysis_group and line.split('\t')[self.columns.index('ANALYSIS_GROUP')] \
                    == analysis_group:
                        data[line.split('\t')[self.columns.index('SAMPLE_NAME')]].\
                        add(line.split('\t')[self.columns.index('RUN_ID')])
                    elif analysis_group and not line.split('\t')\
                    [self.columns.index('ANALYSIS_GROUP')] == analysis_group:
                        continue
                    else:
                        data[line.split('\t')[self.columns.index('SAMPLE_NAME')]]\
                        .add(line.split('\t')[self.columns.index('RUN_ID')])

        return data

