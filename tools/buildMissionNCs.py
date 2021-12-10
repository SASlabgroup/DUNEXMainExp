## buildMission20NC.py - Build Mission 20 Cleaned netCDF
# Import statements
import numpy as np
import sys
import datetime
import logging

# Import DUNEX Tools
sys.path.append('..')
from tools import buildBasicMissionNC
from tools import cleanDUNEXMissions
from tools import missionView
from tools import makeDriftMovie

def main():
    '''
    @edwinrainville

    Description: Build and clean the basic mission datasets.
    TODO: Finish writing description
    '''
    # Set up logger
    logging.basicConfig(filename='../microSWIFT_data/cleanedDataset/build_dataset.log', filemode='w', encoding='utf-8')

    # Define number of missions
    mission_start = 7
    mission_end = 81
    # number_of_missions = 81
    # number_of_missions = 5

    # start timer for entire dataset build
    start_build = datetime.datetime.utcnow()
    
    # Loop through each mission and build netCDF files
    for mission_num in np.arange(mission_start, mission_end+1):
        # Build a netCDF for the mission that has only raw time values read in 
        logging.info('building mission {}'.format(mission_num))
        tic = datetime.datetime.utcnow()
        mission_nc_path = buildBasicMissionNC.main(mission_num=mission_num)
        toc = datetime.datetime.utcnow()
        logging.info('Time to build mission {0} was {1}'.format(mission_num, toc-tic))

        # Clean the dataset
        logging.info('cleaning mission {}'.format(mission_num))
        tic = datetime.datetime.utcnow()
        cleanDUNEXMissions.main(mission_num=mission_num)
        toc = datetime.datetime.utcnow()
        logging.info('Time to clean mission {0} was {1}'.format(mission_num, toc-tic))

        # Plot the mission dataset to view
        logging.info('plotting mission {} dataset'.format(mission_num))
        tic = datetime.datetime.utcnow()
        missionView.main(mission_num=mission_num)
        toc = datetime.datetime.utcnow()
        logging.info('Time to plot mission {0} was {1}'.format(mission_num, toc-tic))

        # # Make Drift movie
        # logging.info('making drfit movie for mission {} dataset'.format(mission_num))
        # tic = datetime.datetime.utcnow()
        # makeDriftMovie.main(mission_num=mission_num)
        # toc = datetime.datetime.utcnow()
        # logging.info('Time to make drift movie for mission {0} was {1}'.format(mission_num, toc-tic))

    # Finish timing the build
    end_build = datetime.datetime.utcnow()
    logging.info('Dataset was built successfully')
    logging.info('Total time to build Dataset = {}'.format(end_build - start_build))


if __name__=='__main__':
    main()
