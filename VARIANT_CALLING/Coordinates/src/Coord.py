'''
Created on 24 Apr 2017

@author: ernesto
'''
import logging

module_logger = logging.getLogger("mainApp.Coord")


class Coord(object):
    '''
    Class to represent a location in a reference
    '''

    def __init__(self, id, start, end):
        '''
        Constructor

         Class variables
        ---------------
        id : str, Required
            ID of the entity for this Coord object (i.e. chr1)
        start : int, Required
            start position of this coordinate in the entity (i.e. 10856277)
        end : int, Required
            end position of this coordinate in the entity (i.e. 10856577)
        '''

        self.id = id
        self.start = start
        self.end = end

    def make_windows(self, step=None):
        '''
        This method will make windows from self.start to self.end (see below).
        The windows created will be end exclusive

        Parameters
        ----------
        step : int , Optional
                step size of each of the windows

        Returns
        ------
        A list containing the different Coord objects, one per window
        '''

        logger = logging.getLogger("mainApp.Coord.make_windows")

        logger.info("Creating windows for id %s and genomic interval "
                    "%d and %d" % (self.id, self.start, self.end))

        ranges = [(n, min(n+step, self.end)) for n in range(self.start, self.end, step)]

        coordlist = []

        for i, item in enumerate(ranges):
            coordlist.append(Coord(id=self.id, start=item[0], end=item[1]-1))

        return coordlist

    def __str__(self):
        sb = []
        for key in self.__dict__:
            sb.append("{key}='{value}'".format(key=key, value=self.__dict__[key]))

        return ', '.join(sb)

    def __repr__(self):
        return self.__str__()
