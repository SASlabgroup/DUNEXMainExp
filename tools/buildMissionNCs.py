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
    logging.basicConfig(filename='../microSWIFT_data/cleanedDataset/build_dataset.log', filemode='w', encoding='utf-8', level='INFO')

    # List of Mission to build for whole dataset
    # skips 6, 47, 53, 55, 57, 64, 65 since they were removed from the dataset
    mission_list = [1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,48,50,51,52,54,56,58,59,60,61,62,63,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81]
    # mission_list = [51]

    # start timer for entire dataset build
    start_build = datetime.datetime.utcnow()
    
    # Loop through each mission and build netCDF files
    for mission_num in mission_list:
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

        # Plot the mission dataset to view after cleaning
        logging.info('plotting mission {} dataset'.format(mission_num))
        tic = datetime.datetime.utcnow()
        filename = 'mission_{}_cleaned.pdf'.format(mission_num)
        missionView.main(mission_num=mission_num, filename=filename)
        toc = datetime.datetime.utcnow()
        logging.info('Time to plot mission {0} was {1}'.format(mission_num, toc-tic)) 
        
       # Make Drift movie after cleaning
        logging.info('making drfit movie for mission {} dataset'.format(mission_num))
        tic = datetime.datetime.utcnow()
        filename = 'mission_{}_drift.gif'.format(mission_num)
        makeDriftMovie.main(mission_num=mission_num, filename=filename)
        toc = datetime.datetime.utcnow()
        logging.info('Time to make drift movie for mission {0} was {1}'.format(mission_num, toc-tic))

    # Finish timing the build
    end_build = datetime.datetime.utcnow()
    logging.info('Dataset was built successfully')
    logging.info('Total time to build Dataset = {}'.format(end_build - start_build))


if __name__=='__main__':
    main()
