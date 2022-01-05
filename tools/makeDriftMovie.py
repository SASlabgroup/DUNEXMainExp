# Make Drift Movie from Mission netCDF
from matplotlib.markers import MarkerStyle
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cftime
import matplotlib.animation as animation

def main(mission_num=None, filename=None):
    '''
    @edwinrainville

    Description: This function will be used to plot a movie of the drift tracks from the DUNEX experiment.
    '''
    # User Inputs    
    if mission_num == None:
        mission_num = int(input('Enter mission number: '))
    
    # define the mission file name
    mission_nc_path = '../microSWIFT_Data/cleanedDataset/mission_{}.nc'.format(mission_num)

    # Create dataset object from the netCDF path
    mission_dataset = nc.Dataset(mission_nc_path, mode='r')

    # Define time from netCDF
    time = cftime.num2pydate(mission_dataset['time'][:], calendar=mission_dataset['time'].calendar, units=mission_dataset['time'].units)

    # Add the FRF Bathymetry to the map 
    # Data from September 28th, 2021
    bathy_url = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/geomorphology/DEMs/surveyDEM/data/FRF_geomorphology_DEMs_surveyDEM_20210928.nc'
    bathy_dataset = nc.Dataset(bathy_url)
    # Create grid from coordinates
    xFRF_grid, yFRF_grid = np.meshgrid(bathy_dataset['xFRF'][:],bathy_dataset['yFRF'][:])
    bathy = bathy_dataset['elevation'][0,:,:]

    # Create map of drift tracks for the microSWIFT
    fig = plt.figure()
    ax = plt.subplot(1,1,1)

    def init_tracks():
        ax.clear()
        ax.set_xlabel('FRF X Location [meters]')
        ax.set_ylabel('FRF Y Location [meters]')
        ax.contourf(xFRF_grid, yFRF_grid, bathy, cmap='gray')

    def update_tracks(index):
        ax.clear()
        ax.set_xlabel('FRF X Location [meters]')
        ax.set_ylabel('FRF Y Location [meters]')
        ax.contourf(xFRF_grid, yFRF_grid, bathy, cmap='gray')

         # Get list of all microSWIFTs on the mission
        microSWIFTs_on_mission = list(mission_dataset.groups.keys())
        for microSWIFT in microSWIFTs_on_mission:
            ax.plot(mission_dataset[microSWIFT]['xFRF'][0], mission_dataset[microSWIFT]['yFRF'][0], color='r', marker='o')
            ax.plot(mission_dataset[microSWIFT]['xFRF'][:index], mission_dataset[microSWIFT]['yFRF'][:index], color='k')
            ax.plot(mission_dataset[microSWIFT]['xFRF'][index], mission_dataset[microSWIFT]['yFRF'][index], color='g', marker='o')
        
        # Add title with time
        ax.set_title('Mission {0} - {1}'.format(mission_num, index))

    # Creating the Animation object
    data_skip = 200
    anim = animation.FuncAnimation(fig, update_tracks, frames=np.arange(0, len(time), data_skip), init_func=init_tracks)
    if filename == None:
        filename = 'mission_{}_drift.gif'.format(mission_num)
    anim.save('../microSWIFT_data/cleanedDataset/Figures/{}'.format(filename), dpi=150, fps=10, writer='pillow')

    # Close dataset and figures
    mission_dataset.close()
    plt.close()
    
if __name__ == '__main__':
    main()