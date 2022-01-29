## Composite Track Map 
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import glob

def main():
    '''
    @edwinrainville

    Description: This function builds a composite map figure of all the tracks of microSWIFTs over the entire expierment to see where we sampled.
    
    '''
    # Get a list of all filenames of missions
    mission_list = glob.glob('../microSWIFT_data/cleanedDataset/mission*.nc')

    # Initialize the plot
    fig_composite_track_map, ax = plt.subplots()
    ax.set_xlabel('FRF X Location [meters]')  
    ax.set_ylabel('FRF Y Location [meters]')
    ax.set_xlim(-100, 900)
    ax.set_ylim(-2000, 2500)

    # For each mission, plot the tracks of each microSWIFT
    for mission in mission_list:
        mission_dataset = nc.Dataset(mission, mode='r')

        # Get list of all microSWIFTs on the mission
        microSWIFTs_on_mission = list(mission_dataset.groups.keys())

         # Loop through all microSWIFTs on mission to build the mission View
        for microSWIFT in microSWIFTs_on_mission:
             # Plot the microSWIFT drift track on bathymetry
            x = mission_dataset[microSWIFT]['xFRF'][:]
            y = mission_dataset[microSWIFT]['yFRF'][:]
            ax.scatter(x, y, marker='.', s=1, linewidths=0)
    
    # Plot the location of our model domain and FRF bathymetry
    ax.plot([0, 800], [-1000, -1000], color='k', label='Current Model Domain')
    ax.plot([0, 800], [2000, 2000], color='k')
    ax.plot([0, 0], [-1000, 2000], color='k')
    ax.plot([800, 800], [-1000, 2000], color='k')
    

    # Plot location of FRF bathymetry
    ax.plot([0, 800], [-100, -100], color='k', linestyle='dashed', label='FRF bathymetry')
    ax.plot([0, 800], [1100, 1100], color='k', linestyle='dashed')
    ax.plot([0, 0], [-100, 1100], color='k', linestyle='dashed')
    ax.plot([800, 800], [-100, 1100], color='k', linestyle='dashed')
    ax.legend()

    # Save the figure 
    plt.savefig('../writing/Figures/CompositeTrackMap.png')
       

if __name__=='__main__':
    main()