## buildMission20NC.py - Build Mission 20 Cleaned netCDF
# Import statements
import numpy as np
import sys
import datetime
import logging

# Import DUNEX Tools
sys.path.append('..')
from tools import buildBasicMissionNC

def main():
    '''
    @edwinrainville

    Description: Build and clean the basic mission datasets.
    TODO: Finish writing description
    '''
    # Set up logger
    logging.basicConfig(filename='../microSWIFT_data/cleanedDataset/build_dataset.log', encoding='utf-8', level=logging.DEBUG)

    # Define number of missions
    number_of_missions = 81
    
    # Loop through each mission and build netCDF files
    for mission_num in np.arange(1, number_of_missions):
        # Build a netCDF for the mission that has only raw time values read in 
        logging.info('building mission {}'.format(mission_num))
        tic = datetime.datetime.utcnow()
        mission_nc_path = buildBasicMissionNC.main(mission_num=mission_num)
        toc = datetime.datetime.utcnow()
        logging.info('Time to build mission {0} was {1}'.format(mission_num, toc-tic))

if __name__=='__main__':
    main()
