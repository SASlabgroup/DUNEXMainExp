# Make Drift Movie from Mission netCDF
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cftime
import matplotlib.animation as animation

def main(mission_num=None):
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

    # Get list of all microSWIFTs on the mission
    microSWIFTs_on_mission = list(mission_dataset.groups.keys())
    microSWIFT =microSWIFTs_on_mission[0]

    # Add the FRF Bathymetry to the map 
    # Data from September 28th, 2021
    bathy_url = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/geomorphology/DEMs/surveyDEM/data/FRF_geomorphology_DEMs_surveyDEM_20210928.nc'
    bathy_dataset = nc.Dataset(bathy_url)
    # Create grid from coordinates
    xFRF_grid, yFRF_grid = np.meshgrid(bathy_dataset['xFRF'][:],bathy_dataset['yFRF'][:])
    bathy = bathy_dataset['elevation'][0,:,:]

    # Create map of drift tracks for the microSWIFT
    fig = plt.figure()
    # fig.set_size_inches(8.5, 11)
    ax = plt.subplot(1,1,1)

    def init_tracks():
        ax.clear()
        ax.set_xlabel('FRF X Location [meters]')
        ax.set_ylabel('FRF Y Location [meters]')
        ax.contourf(xFRF_grid, yFRF_grid, bathy, cmap='gray')
        
    def update_tracks(index):
        ax.scatter(mission_dataset[microSWIFT]['xFRF'][index], mission_dataset[microSWIFT]['yFRF'][index])

    # Creating the Animation object
    anim = animation.FuncAnimation(fig, update_tracks, frames=np.arange(0, len(time)), init_func=init_tracks, interval=20)
    anim.save('drifttrack.mp4', dpi=150, fps=30, writer='ffmpeg')

    # Close dataset and figures
    mission_dataset.close()
    plt.close()
    
if __name__ == '__main__':
    main()