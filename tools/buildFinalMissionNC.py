"""
Build the final DUNEX microSWIFT dataset. 
"""
import datetime
import os

import netCDF4 as nc
import numpy as np

def main():
    """

    """
    
    buildFinalNC(mission_num=1)

def buildFinalNC(mission_num):
    """
    Function that builds a netCDF file from the varaibles processed in
    each mission. These netCDF files are compliant with the CF-1.6 
    conventions. Each microSWIFT track will be treated as a trajectory
    with a trajectory id for the microSWIFT number.
    """
    try:
        mission_nc_path = '../microSWIFT_data/cleanedDataset/' \
                          'mission_{}.nc'.format(mission_num)
        cleaned_dataset = nc.Dataset(mission_nc_path, mode='r')
    except:
        raise FileNotFoundError("Can't find this mission file.")

    
    # Create new file for the final dataset
    final_data_fname = '../microSWIFT_data/final_dataset/' \
                       'mission_{}.nc'.format(mission_num)

    # Open netCDF dataset
    rootgrp = nc.Dataset(final_data_fname, 'w',
                         clobber=True, 
                         format='NETCDF4_CLASSIC')

    # Set overall file information
    rootgrp.Conventions = 'CF-1.6'
    rootgrp.title = 'DUNEX microSWIFT drifter - Mission {}'.format(mission_num)
    rootgrp.institution = 'University of Washington - Applied Physics Lab'
    rootgrp.source = 'Observations from microSWIFT drifters deployed ' \
                     'during DUNEX'
    rootgrp.history = str(datetime.datetime.utcnow()) + ' Python'
    rootgrp.references = 'https://github.com/SASlabgroup/microSWIFT and ' \
                         'https://github.com/SASlabgroup/DUNEXMainExp'
    rootgrp.comment = get_cleaning_notes(mission_num)

    # Create the dimensions for the dataset
    mission_time_dim = rootgrp.createDimension('time', cleaned_dataset['time'][:].shape[0])
    mission_time_nc = rootgrp.createVariable('time', np.float64, ('time',))
    mission_time_nc.units = cleaned_dataset['time'].units
    mission_time_nc.calendar = cleaned_dataset['time'].calendar
    mission_time_nc.standard_name = 'time'
    mission_time_nc.long_name = 'time'
    mission_time_nc[:] = cleaned_dataset['time'][:]

    microSWIFTs_on_mission = list(cleaned_dataset.groups.keys())
    microSWIFT_nums = [int(microSWIFT[11:]) for microSWIFT in microSWIFTs_on_mission]
    trajectory_dim = rootgrp.createDimension('trajectory', len(microSWIFTs_on_mission))
    trajectory_nc = rootgrp.createVariable('trajectory', np.intc, ('trajectory',))
    trajectory_nc.cf_role = 'trajectory_id'
    trajectory_nc.long_name = 'trajectory name'
    trajectory_nc[:] = microSWIFT_nums

    # Add sampling frequency information to the netcdf file
    single_val_dim = rootgrp.createDimension('single_value', 1)
    gps_freq = rootgrp.createVariable('gps_freq', np.intc, ('single_value'))
    gps_freq.long_name = 'GPS module sampling frequency'
    gps_freq.units = 's-1'
    gps_freq[:] = cleaned_dataset['gps_freq'][:]
    
    imu_freq = rootgrp.createVariable('imu_freq', np.intc, ('single_value'))
    imu_freq.long_name = 'IMU module sampling frequency'
    imu_freq.units = 's-1'
    imu_freq[:] = cleaned_dataset['imu_freq'][:]

    # Body Accelerations
    accel_x_body_nc = rootgrp.createVariable('accel_x_body',
                                             np.float64,
                                             ('trajectory', 'time'))
    accel_x_body_nc.long_name = 'X-axis acceleration in bouy reference frame'
    accel_x_body_nc.description = 'Accleration across the minor axis of the ' \
                                  'buoy in the reference frame of the bouy.'
    accel_x_body_nc.units = 'm s-2'
    save_var_in_final_style('accel_x_body',
                            accel_x_body_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)         

    rootgrp.close()
    

    return 

def get_cleaning_notes(mission_num):
    """
    Get the cleaning notes for each mission.
    """
    comment = 'yay this is a comment'

    return comment

def save_var_in_final_style(var_name,
                            nc_var,
                            microSWIFTs_on_mission,
                            cleaned_dataset):
    """
    Save the given variable in the new nc format from old file.

    Parameters
    ----------
    var_name : _type_
        _description_
    """
    for n in range(len(microSWIFTs_on_mission)):
        nc_var[n,:] = cleaned_dataset[microSWIFTs_on_mission[n]][var_name][:]      

if __name__ == '__main__':
    main()