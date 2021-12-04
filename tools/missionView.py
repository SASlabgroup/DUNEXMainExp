## missionView.py
# Import Packages
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cftime
import matplotlib.backends.backend_pdf

def main(mission_nc_path=None, mission_num=None):
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
    # User Inputs
    if mission_nc_path == None:
        mission_nc_path = input('Enter path to mission netCDF:')
    
    if mission_num == None:
        mission_num = int(input('Enter mission number:'))

    # Create dataset object from the netCDF path
    mission_dataset = nc.Dataset(mission_nc_path, mode='r')

    # Define time from netCDF
    time = cftime.num2pydate(mission_dataset['time'][:], calendar=mission_dataset['time'].calendar, units=mission_dataset['time'].units)

    # Set up pdf to save figures to for each microSWIFT
    pdf_name = '../buildDataset/missionViews/mission_{}.pdf'.format(mission_num)
    pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_name)

    # Get list of all microSWIFTs on the mission
    microSWIFTs_on_mission = list(mission_dataset.groups.keys())

    # Loop through all microSWIFTs on mission to build the mission View
    for microSWIFT_num in np.arange(len(microSWIFTs_on_mission)):

        # Create map of drift tracks for the microSWIFT
        fig = plt.figure()
        fig.set_size_inches(8.5, 11)
        ax1 = fig.add_subplot(1,2,1)
        ax1.set_xlabel('FRF X Location [meters]')  
        ax1.set_ylabel('FRF Y Location [meters]')

        # Add the FRF Bathymetry to the map 
        # Data from September 28th, 2021
        bathy_url = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/geomorphology/DEMs/surveyDEM/data/FRF_geomorphology_DEMs_surveyDEM_20210928.nc'
        bathy_dataset = nc.Dataset(bathy_url)
        # Create grid from coordinates
        xFRF_grid, yFRF_grid = np.meshgrid(bathy_dataset['xFRF'][:],bathy_dataset['yFRF'][:])
        bathy = bathy_dataset['elevation'][0,:,:]
        ax1.contourf(xFRF_grid, yFRF_grid, bathy, cmap='gray')

        # Plot the microSWIFT drift track on bathymetry
        x = mission_dataset[microSWIFTs_on_mission[microSWIFT_num]]['xFRF'][:]
        y = mission_dataset[microSWIFTs_on_mission[microSWIFT_num]]['yFRF'][:]
        ax1.scatter(x, y, color='g')
        ax1.set_xlim([np.min(x)-100, np.max(x)+100])
        ax1.set_ylim([np.min(y)-100, np.max(y)+100])
        ax1.set_title('Drift Track')

        # Plot the accelerations
        ax2 = fig.add_subplot(3,2,2)
        ax2.plot(time, mission_dataset[microSWIFTs_on_mission[microSWIFT_num]]['accel_x'][:], color='g', label='X')
        ax2.plot(time, mission_dataset[microSWIFTs_on_mission[microSWIFT_num]]['accel_y'][:], color='b', label='Y')
        ax2.plot(time, mission_dataset[microSWIFTs_on_mission[microSWIFT_num]]['accel_z'][:], color='k', label='Z')
        ax2.set_xlabel('Time')
        ax2.set_ylabel('Acceleration [m/s^2]')
        ax2.legend(bbox_to_anchor=(0,1.04,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=3)

        # Plot the Gyroscope
        ax3 = fig.add_subplot(3,2,4)
        ax3.plot(time, mission_dataset[microSWIFTs_on_mission[microSWIFT_num]]['gyro_x'][:], color='g', label='X')
        ax3.plot(time, mission_dataset[microSWIFTs_on_mission[microSWIFT_num]]['gyro_y'][:], color='b', label='Y')
        ax3.plot(time, mission_dataset[microSWIFTs_on_mission[microSWIFT_num]]['gyro_z'][:], color='k', label='Z')
        ax3.set_xlabel('Time')
        ax3.set_ylabel('Rotations [degrees/sec]')
        ax3.set_title('Gyroscope')

        # Plot the Magnetometer
        ax4 = fig.add_subplot(3,2,6)
        ax4.plot(time, mission_dataset[microSWIFTs_on_mission[microSWIFT_num]]['mag_x'][:], color='g', label='X')
        ax4.plot(time, mission_dataset[microSWIFTs_on_mission[microSWIFT_num]]['mag_y'][:], color='b', label='Y')
        ax4.plot(time, mission_dataset[microSWIFTs_on_mission[microSWIFT_num]]['mag_z'][:], color='k', label='Z')
        ax4.set_xlabel('Time')
        ax4.set_ylabel('Magnetometer [uTeslas]')
        ax4.set_title('Magnetometer')

        # Figure Properties 
        plt.tight_layout()
        ax1.set_aspect('equal')
        pdf.savefig(fig)

    # Close the dataset
    mission_dataset.close()
    bathy_dataset.close()
    pdf.close()

if __name__=='__main__':
    main()