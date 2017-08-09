from Coord import *
import logging

if __name__ == '__main__':
    
    logger = logging.getLogger("mainApp")
    logger.setLevel(logging.INFO)

    # create the logging file handler
    fh = logging.FileHandler("create_windows.log")

    # add handler to logger object
    logger.addHandler(fh)
 
    logger.info("Program started")

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)

    c=Coord(id='chr1',start=0,end=1000)
    l=c.make_windows(step=100)

    logger.info("Done!")

