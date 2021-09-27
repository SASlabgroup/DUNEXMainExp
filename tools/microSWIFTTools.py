# Tools for analyzing data from microSWIFTs
# Import modules
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import netCDF4 as nc


def missionMap(mission_num, mission_dir_path, mission_nc_path):
    '''
    @edwinrainville

    '''
    # Load in netCDF file as a dataset
    mission_dataset = nc.Dataset(mission_nc_path, mode='r')
    
    # Get list of all microSWIFTs on the mission
    microSWIFTs_on_mission = list(mission_dataset.groups.keys())

    # Create map of all drift tracks during mission
    fig, ax = plt.subplots()
    ax.set_xlabel('Longitude')  
    ax.set_ylabel('Latitude')
    fig.autofmt_xdate()
    for microSWIFT in microSWIFTs_on_mission:
        ax.scatter(mission_dataset[microSWIFT]['GPS']['lon'][:], mission_dataset[microSWIFT]['GPS']['lat'][:], c=mission_dataset[microSWIFT]['GPS']['time'][:])

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
    
