# Clean each DUNEX mission
import numpy as np
import logging
import netCDF4 as nc
import matplotlib.pyplot as plt
import cftime

def computeBeachLocation():
    '''
    This function creates a time series of furthest offshore point of dry beach from bathymetry and 

    Returns:
    time series of cross shore extent of beach 

    '''
    # load in bathymetry 
    url = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/geomorphology/DEMs/surveyDEM/data/FRF_geomorphology_DEMs_surveyDEM_20211021.nc'
    bathy = nc.Dataset(url)
    xFRF = bathy['xFRF'][:]
    yFRF = bathy['yFRF'][:]
    xFRF_grid, yFRF_grid = np.meshgrid(xFRF, yFRF)
    elevation = bathy['elevation'][0,:,:]

    # Load in water level data
    # TODO: replace these values with true waterLavel Data instead of predicted once it is posted
    url = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/oceanography/waterlevel/eopNoaaTide/2021/FRF-ocean_waterlevel_eopNoaaTide_202110.nc'
    waterLevel_dataset = nc.Dataset(url)
    time = cftime.num2pydate(waterLevel_dataset['time'][:], calendar=waterLevel_dataset['time'].calendar, units=waterLevel_dataset['time'].units)
    # waterLevel_series = waterLevel_dataset['waterLevel'][:]
    waterLevel_series = waterLevel_dataset['predictedWaterLevel'][:]

    # Define empty list for max x values
    beach_x_series = []
    
    # define cross shore buffer
    buffer = 5 # units are meters

    for waterLevel in waterLevel_series:

        # Compute depth 
        depth = elevation + waterLevel
        depth[depth < 0] = 0
        depth[depth > 0 ] = 1

        # define empty list to fill with values for beach x location
        max_x_row = []
        for n in np.arange(len(yFRF)):
            max_x_row.append(xFRF[np.where(depth[n,:] == 0)[0][0]])
        max_x = np.max(max_x_row)
        max_x += buffer

        # Append to time series
        beach_x_series.append(max_x)

    # fig, ax = plt.subplots()
    # ax.plot(beach_x_series)
    # plt.savefig('max_x.png')

    # return the time series
    return beach_x_series, time

def individualMissionMicroSWIFTMask(mission_num = None):
    '''
    Returns the indivudal mask for each microSWIFT in each mission which are developed from visual inspection - each mask has notes as to how it was developed
    TODO: Make spreadsheet of all individual mask slices
    '''
    # Define mission Number
    if mission_num == None:
        mission_num = int(input('Enter Mission Number: '))

    # Function to initialize mask dictionary
    def init_mask_dict(mission_num=mission_num):
        # define mission netCF
        mission_dataset = nc.Dataset('../microSWIFT_data/cleanedDataset/mission_{}.nc'.format(mission_num))

        # define list of microSWIFTs on mission 
        microSWIFTs_on_mission = list(mission_dataset.groups.keys())

        # Initialize mask 
        mission_masks = {microSWIFT: None for microSWIFT in microSWIFTs_on_mission}

        # Return the initialized mask
        return mission_masks
    
    # Mission 1 Individual mask 
    if mission_num == 1:
        # Intialize mask 
        mission_masks = init_mask_dict(mission_num=mission_num)

        # Append Mask with indices to mask

    # Mission 2 Individual mask 
    if mission_num == 2:
        # Intialize mask 
        mission_masks = init_mask_dict(mission_num=mission_num)

        # Append Mask with indices to mask

    # Mission 3 Individual mask 
    if mission_num == 3:
        # Intialize mask 
        mission_masks = init_mask_dict(mission_num=mission_num)

        # Append Mask with indices to mask

    # Mission 4 Individual mask 
    if mission_num == 4:
        # Intialize mask 
        mission_masks = init_mask_dict(mission_num=mission_num)

        # Append Mask with indices to mask

    # Mission 5 Individual mask 
    if mission_num == 5:
        # Intialize mask 
        mission_masks = init_mask_dict(mission_num=mission_num)

        # Append Mask with indices to mask

    # Mission 6 Individual mask 
    if mission_num == 6:
        # Intialize mask 
        mission_masks = init_mask_dict(mission_num=mission_num)

        # Append Mask with indices to mask

    # Mission 7 Individual mask 
    if mission_num == 7:
        # Intialize mask 
        mission_masks = init_mask_dict(mission_num=mission_num)

        # Append Mask with indices to mask

    # Mission 8 Individual mask 
    if mission_num == 8:
        # Intialize mask 
        mission_masks = init_mask_dict(mission_num=mission_num)

        # Append Mask with indices to mask


    # Mission 9 Individual mask 
    if mission_num == 9:
        # Intialize mask 
        mission_masks = init_mask_dict(mission_num=mission_num)

        # Append Mask with indices to mask

    # Mission 10 Individual mask 
    if mission_num == 10:
        # Intialize mask 
        mission_masks = init_mask_dict(mission_num=mission_num)

        # Append Mask with indices to mask

     # Mission 20 Individual mask 
    if mission_num == 20:
        # Intialize mask 
        mission_masks = init_mask_dict(mission_num=mission_num)

        # microSWIFT_31
        # description: There is a section of all imu data at the end of the time series that values are not fluctuating signifying that it is on the beach but not masked by the 
        mission_masks['microSWIFT_31'] = slice(9000, -1)

        # microSWIFT_21
        # description: There is a very small section of data at the very end of the time series that it most likely is being masked from the beach but maybe moved back across the beach mask briefly 
        mission_masks['microSWIFT_21'] = slice(10000, -1)
        
    # Return the mission mask dictionary
    return mission_masks

def main(mission_num=None):
    '''
    @edwinrainville

    Description: Clean each DUNEX dataset by masking invalid points

    '''
    # Set up logger
    logging.basicConfig(filename='../microSWIFT_data/cleanedDataset/build_dataset.log', encoding='utf-8', level=logging.DEBUG)

     # Compute beach 
    beach_x_series, beach_x_time = computeBeachLocation()

    # Load in each mission netCDF and apply masks 
    data_dir = '../microSWIFT_data/cleanedDataset/'

    # Define mission Number
    if mission_num == None:
        mission_num = int(input('Enter Mission Number: '))

    # Convert mission number to a path to the file
    mission_nc = 'mission_{}.nc'.format(mission_num)
    mission_nc_path = data_dir + mission_nc

    # Load in mission datset as a netCDF object
    mission_dataset = nc.Dataset(mission_nc_path, mode='a')

    # From start time find the correct beach_x value 
    mission_start = cftime.num2pydate(mission_dataset['time'][0], calendar=mission_dataset['time'].calendar, units=mission_dataset['time'].units)
    start_time_closest_index = int(np.argmin(np.abs(beach_x_time - mission_start)))
    beach_x_during_mission = beach_x_series[start_time_closest_index]

    # For each microSWIFT on mission
    microSWIFTs_on_mission = list(mission_dataset.groups.keys())
    for microSWIFT in microSWIFTs_on_mission:
        # Get values for cross shore location of microSWIFT
        xFRF = mission_dataset[microSWIFT]['xFRF'][:]

        # Find first index in the xFRF data that is less than the beach location
        beach_mask_ind = np.where(xFRF <= beach_x_during_mission)

        # Mask all indices where the microSWIFT is on the beach
        microSWIFT_variables = list(mission_dataset[microSWIFT].variables.keys())
        for variable in microSWIFT_variables:
            mission_dataset[microSWIFT][variable][beach_mask_ind] = np.ma.masked

    # Mask each individual microSWIFT on each mission
    mission_mask_dict = individualMissionMicroSWIFTMask(mission_num=mission_num)
    for microSWIFT in list(mission_mask_dict.keys()):
        microSWIFT_mask = mission_mask_dict[microSWIFT]
        if microSWIFT_mask != None:
            # Mask all indices where the microSWIFT is on the beach
            microSWIFT_variables = list(mission_dataset[microSWIFT].variables.keys())
            for variable in microSWIFT_variables:
                mission_dataset[microSWIFT][variable][microSWIFT_mask] = np.ma.masked
    
    # Close the dataset
    mission_dataset.close()


if __name__ == '__main__':
    main()