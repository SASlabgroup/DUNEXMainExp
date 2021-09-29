# Tools for analyzing data from microSWIFTs
# Import modules
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import cftime
import datetime


def missionMap(mission_num, mission_dir_path, mission_nc_path):
    '''
    @edwinrainville

    '''
    # Load in netCDF file as a dataset
    mission_dataset = nc.Dataset(mission_nc_path, mode='r')
    
    # Get list of all microSWIFTs on the mission
    microSWIFTs_on_mission = list(mission_dataset.groups.keys())

    # Create map of all drift tracks during mission
    fig, ax = plt.subplots(figsize=(8,6))
    ax.set_xlabel('Longitude')  
    ax.set_ylabel('Latitude')
    fig.autofmt_xdate()

    # Sort time labels 
    min_time_label = mission_dataset[microSWIFTs_on_mission[0]]['GPS']['time'][0]
    max_time_label = mission_dataset[microSWIFTs_on_mission[0]]['GPS']['time'][-1]

    for microSWIFT in microSWIFTs_on_mission:
        # Plot the lat-lon map with color coded points in time
        map = ax.scatter(mission_dataset[microSWIFT]['GPS']['lon'][:], mission_dataset[microSWIFT]['GPS']['lat'][:], c=mission_dataset[microSWIFT]['GPS']['time'][:])
        # reset time labels for colormap 
        if mission_dataset[microSWIFT]['GPS']['time'][0] < min_time_label:
            min_time_label = mission_dataset[microSWIFT]['GPS']['time'][0]
        if mission_dataset[microSWIFT]['GPS']['time'][-1] > max_time_label:
            max_time_label = mission_dataset[microSWIFT]['GPS']['time'][-1]
    
    # Set colorbar and figure properties
    cbar = fig.colorbar(map, ax=ax, ticks=[min_time_label, max_time_label])
    map.set_clim([min_time_label, max_time_label])
    cbar.ax.set_xlabel('Time [UTC]')
    time_labels = cftime.num2pydate([min_time_label, max_time_label],units=mission_dataset[microSWIFT]['GPS']['time'].units, calendar=mission_dataset[microSWIFT]['GPS']['time'].calendar)
    time_labels = [time_labels[0].strftime('%Y-%m-%d %H:%M'), time_labels[1].strftime('%Y-%m-%d %H:%M')]
    cbar.ax.set_yticklabels(time_labels, rotation=0, va='center')
    plt.tight_layout()
    figure_path = mission_dir_path + 'mission_{}_map.png'.format(mission_num)
    plt.savefig(figure_path)
    return figure_path

def missionAccels():
    '''
    @edwinrainville
    '''
    return 'This function plots a time series of all accelerations '

def computeWaveProperties():
    '''
    @edwinrainville

    Description: This function computes all bulk wave properties about the field from the drifter tracks during a microSWIFT 
    mission.

    '''

def reconstructWaveField(mission_dir_path):
    '''
    @edwinrainville

    Description: This funnction uses a compressive sensing algorithm to reconstruct a field of waves during a 
    microSWIFT mission. 
    '''
    
def FRFtransform(lat,lon):
    '''
    @edwinrainville, Originally written by J. Thomson, 1/2011

    Description: function to convert from lat & lon (decimal degrees, negative longitude) to FRF x,y (meters)
    '''

    # Define offsets
    latoffset = 36.178039
    lonoffset = -75.749672

    # Define constants
    rotation = 19 #rotation in degress CCW from True north
    radius = 6371*np.cos(np.radians(36)) # what is 36?

    

    x = ''
    y = ''
    return x, y