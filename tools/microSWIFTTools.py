# Tools for analyzing data from microSWIFTs
# Import modules
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import netCDF4 as nc


def mission_map(mission_num, mission_dir_path, mission_nc_path):
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

def mission_accels():
    '''
    @edwinrainville
    '''
    return 'This function plots a time series of all accelerations '
    
