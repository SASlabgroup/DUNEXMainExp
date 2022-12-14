# Clean each DUNEX mission
import numpy as np
import logging
import netCDF4 as nc
import matplotlib.pyplot as plt
import cftime
import pandas as pd

def computeBeachLocation():
    '''
    This function creates a time series of furthest offshore point of dry beach from bathymetry and 

    Returns:
    time series of cross shore extent of beach 

    '''
    # load in bathymetry 
    # path = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/geomorphology/DEMs/surveyDEM/data/FRF_geomorphology_DEMs_surveyDEM_20211021.nc'
    path = '../microSWIFT_data/FRFdata/FRF_geomorphology_DEMs_surveyDEM_20211021.nc'
    bathy = nc.Dataset(path)
    xFRF = bathy['xFRF'][:]
    yFRF = bathy['yFRF'][:]
    xFRF_grid, yFRF_grid = np.meshgrid(xFRF, yFRF)
    elevation = bathy['elevation'][0,:,:]

    # Load in water level data
    # url = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/oceanography/waterlevel/eopNoaaTide/2021/FRF-ocean_waterlevel_eopNoaaTide_202110.nc'
    path = '../microSWIFT_data/FRFdata/FRF-ocean_waterlevel_eopNoaaTide_202110.nc'
    waterLevel_dataset = nc.Dataset(path)
    time = cftime.num2pydate(waterLevel_dataset['time'][:], calendar=waterLevel_dataset['time'].calendar, units=waterLevel_dataset['time'].units)
    # waterLevel_series = waterLevel_dataset['waterLevel'][:]
    waterLevel_series = waterLevel_dataset['waterLevel'][:]

    # Define empty list for max x values
    beach_x_series = []
    
    # define cross shore buffer
    buffer = 2 # units are meters

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

    # return the time series
    return beach_x_series, time

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

    if mission_num == 51 or mission_num == 52:
        beach_x_during_mission = -100 # for mission 51 and 52 since they have been cleaned very throoughyl with hand written notes

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

    # Read in mask values from spreadsheet and mask all those values
    dunex_notes_filename = '../DUNEXMainExp_notes.xlsx'
    dunex_df = pd.read_excel(dunex_notes_filename, 'microSWIFT Masks')
    dunex_df['Mission Number'] = pd.Series(dunex_df['Mission Number']).fillna(method='ffill')

    # Find all rows in mission 
    mission_masks = dunex_df[dunex_df['Mission Number'] == mission_num]

    # microSWIFT Masks
    microSWIFTs_masked_on_mission = list(dunex_df[dunex_df['Mission Number'] == mission_num]['microSWIFT ID'])
    if type(microSWIFTs_masked_on_mission[0]) == str:
        for microSWIFT in microSWIFTs_masked_on_mission:
            # Get list of microSWIFT Varaibles
            try:
                microSWIFT_variables = list(mission_dataset[microSWIFT].variables.keys())

                # Get the mask values from spreadsheet 
                individual_microSWIFT_mask = mission_masks[mission_masks['microSWIFT ID'] == microSWIFT]

                # Mask from Begininng to start mask end index 
                start_mask_end_index = int(individual_microSWIFT_mask['Start Mask End Index'].item())
                end_mask_start_index = int(individual_microSWIFT_mask['End Mask Start Index'].item())
                slice_list = []
                additional_masked_points = []
                if individual_microSWIFT_mask['Additional Masking indices'].item() != 0:
                    # sort between ranges we want to mask and individual points
                    for val in individual_microSWIFT_mask['Additional Masking indices'].item().split(','):
                        if ':' in val:
                            slice_list.append(slice(int(val.split(':')[0]), int(val.split(':')[1])))
                        else:
                            additional_masked_points.append(int(val))

                    # Convert additional points list to numpy array
                    additional_masked_points = np.array(additional_masked_points)

                # Mask all indices in the mask list
                for variable in microSWIFT_variables:
                    mission_dataset[microSWIFT][variable][:start_mask_end_index] = np.ma.masked
                    mission_dataset[microSWIFT][variable][end_mask_start_index:] = np.ma.masked
                    if len(additional_masked_points) > 0:
                        mission_dataset[microSWIFT][variable][additional_masked_points] = np.ma.masked
                    if len(slice_list) > 0:
                        for inds in slice_list:
                            mission_dataset[microSWIFT][variable][inds] = np.ma.masked
            except: 
                continue

    # Close the dataset
    mission_dataset.close()

if __name__ == '__main__':
    main()