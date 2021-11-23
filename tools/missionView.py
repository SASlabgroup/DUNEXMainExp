## missionView.py
# Import Packages
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cftime

def main(mission_nc_path, mission_num):
    '''
    @edwinrainville

    Description: A quick view of all data in a mission netCDF.

    Parameters:
        mission_nc: file name of mission netCDF file

    Returns:
        Drift tracks: Scatter plots of x and y position of the microSWIFT

        Accelerations:

        Magnetometers: 

        Gyroscope: 

        Click Through each microSWIFT in a mission 
    
    '''
    # Create dataset object from the netCDF path
    mission_dataset = nc.Dataset(mission_nc_path, mode='r')
    
    # Get list of all microSWIFTs on the mission
    microSWIFTs_on_mission = list(mission_dataset.groups.keys())

    # Create map of all drift tracks during mission
    fig, ax = plt.subplots(figsize=(8,6))
    ax.set_xlabel('FRF X Location [meters]')  
    ax.set_ylabel('FRF Y Location [meters]')

    # Add the FRF Bathymetry to the map 
    # Data from September 28th, 2021
    bathy_url = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/geomorphology/DEMs/surveyDEM/data/FRF_geomorphology_DEMs_surveyDEM_20210928.nc'
    bathy_dataset = nc.Dataset(bathy_url)
    # Create grid from coordinates
    xFRF_grid, yFRF_grid = np.meshgrid(bathy_dataset['xFRF'][:],bathy_dataset['yFRF'][:])
    bathy = bathy_dataset['elevation'][0,:,:]
    ax.contourf(xFRF_grid, yFRF_grid, bathy, cmap='gray')

    # # Sort time labels 
    # initial_values_set = False
    # while initial_values_set == False:
    #     for microSWIFT in microSWIFTs_on_mission:
    #         if 'GPS' in list(mission_dataset[microSWIFT].groups.keys()):
    #             if 'time' in list(mission_dataset[microSWIFT]['GPS'].variables):
    #                 # Set initial time labels
    #                 min_time_label = mission_dataset[microSWIFT]['GPS']['time'][0]
    #                 max_time_label = mission_dataset[microSWIFT]['GPS']['time'][-1]

    #                 # Sort x and y locations for map  
    #                 min_x = np.min(mission_dataset[microSWIFT]['GPS']['xFRF'])
    #                 max_x = np.max(mission_dataset[microSWIFT]['GPS']['xFRF'])
    #                 min_y = np.min(mission_dataset[microSWIFT]['GPS']['yFRF'])
    #                 max_y = np.max(mission_dataset[microSWIFT]['GPS']['yFRF'])
    #                 initial_values_set = True
    #             else:
    #                 continue
    #         else:
    #             continue
                
    # for microSWIFT in microSWIFTs_on_mission:
    #     # Compute local coordinates from each lat-lon series
    #     if 'GPS' in list(mission_dataset[microSWIFT].groups.keys()):
    #         if 'lat' and 'lon' in list(mission_dataset[microSWIFT]['GPS'].variables):
    #             # Plot the lat-lon map with color coded points in time
    #             map = ax.scatter(mission_dataset[microSWIFT]['GPS']['xFRF'], mission_dataset[microSWIFT]['GPS']['yFRF'], c=mission_dataset[microSWIFT]['GPS']['time'][:], cmap='plasma')
    #             # reset time labels for colormap 
    #             if mission_dataset[microSWIFT]['GPS']['time'][0] < min_time_label:
    #                 min_time_label = mission_dataset[microSWIFT]['GPS']['time'][0]
    #             if mission_dataset[microSWIFT]['GPS']['time'][-1] > max_time_label:
    #                 max_time_label = mission_dataset[microSWIFT]['GPS']['time'][-1]
    #             # Set max and min position
    #             if np.min(x) < min_x:
    #                 min_x = np.min(x)
    #             if np.max(x) > max_x:
    #                 max_x = np.max(x)
    #             if np.min(y) < min_y:
    #                 min_y = np.min(y)
    #             if np.max(y) > max_y:
    #                 max_y = np.max(y)
                
    #             # Figure Properties
    #             time_labels = cftime.num2pydate([min_time_label, max_time_label],units=mission_dataset[microSWIFT]['GPS']['time'].units, calendar=mission_dataset[microSWIFT]['GPS']['time'].calendar)
            
    # # Set colorbar and figure properties
    # cbar = fig.colorbar(map, ax=ax, ticks=[min_time_label, max_time_label])
    # map.set_clim([min_time_label, max_time_label])
    # ax.set_xlim([min_x, max_x])
    # ax.set_ylim([min_y, max_y])
    # cbar.ax.set_xlabel('Time [UTC]')
    # time_labels = [time_labels[0].strftime('%Y-%m-%d %H:%M'), time_labels[1].strftime('%Y-%m-%d %H:%M')]
    # cbar.ax.set_yticklabels(time_labels, rotation=0, va='center')

    plt.tight_layout()
    ax.set_aspect('equal')
    figure_path = './mission_{}_map.png'.format(mission_num)
    plt.savefig(figure_path)

    # Close the dataset
    mission_dataset.close()
    bathy_dataset.close()

if __name__=='__main__':
    main()