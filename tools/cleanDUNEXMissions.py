# Clean each DUNEX mission
import numpy as np
import logging
import netCDF4 as nc

def computeBeachLocation():
    '''
    This function creates a time series of furthest offshore point of dry beach from bathymetry and 

    Returns:
    time series of cross shore extent of beach 

    '''
    # load in bathymetry 
    url = ''
    bathy = nc.Dataset(url)

    # Load in water level data
    waterLevel = 0  # units are meters

    # define empty list to fill with values for beach x location
    beach_x_series = []

    # Compute depth 
    depth = bathy + waterLevel
    depth[depth < 0] = 0
    depth[depth > 0 ] = 1
    # find last index that is one for each row

    # Find max last index and then the beach xFRF location
    beach_x = 20

    # define cross shore buffer
    buffer = 20 # units are meters

    # Append to time series
    beach_x_series.append(beach_x + buffer)

    # return the time series
    return beach_x_series

def main():
    '''
    @edwinrainville

    Description: Clean each DUNEX dataset by masking invalid points

    '''
    # Set up logger
    logging.basicConfig(filename='../microSWIFT_data/cleanedDataset/build_dataset.log', encoding='utf-8', level=logging.DEBUG)

     # Compute beach 
    beach_x_series = computeBeachLocation()

    # Define number of missions
    number_of_missions = 81

    # Load in each mission netCDF and apply masks 
    data_dir = '../microSWIFT_data/cleanedDataset/'

    # Define mission Number
    mission_num = 20
    mission_nc = 'mission_{}.nc'.format(mission_num)

    # Load in mission datset as a netCDF object
    mission_dataset = nc.Dataset(mission_nc, mode='a')

    # For each microSWIFT on mission
    microSWIFTs_on_mission = list(mission_dataset.groups.keys())
    # for microSWIFT in microSWIFTs_on_mission:
    microSWIFT = microSWIFTs_on_mission[0]

    # Close the dataset
    mission_dataset.close()





if __name__ == '__main__':
    main()