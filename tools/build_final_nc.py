"""
Build the final DUNEX microSWIFT dataset.
"""
import datetime

import pandas as pd
import netCDF4 as nc
import numpy as np

def main():
    """
    Runs the buildFinalNC function to test.
    """
    # missions = [1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
    #             23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,
    #             41,42,43,44,45,46,48,50,51,52,54,56,58,59,60,61,62,63,
    #             66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81]
    missions = [78]
    for mission in missions:
        build_final_nc(mission_num=mission)

def build_final_nc(mission_num):
    """
    Function that builds a netCDF file from the variables processed in
    each mission. These netCDF files are compliant with the CF-1.6
    conventions. Each microSWIFT track will be treated as a trajectory
    with a trajectory id for the microSWIFT number.
    """
    try:
        mission_nc_path = '../microSWIFT_data/cleanedDataset/' \
                          f'mission_{mission_num}.nc'
        cleaned_dataset = nc.Dataset(mission_nc_path, mode='r')
    except Exception as exc:
        raise FileNotFoundError("Can't find this mission file.") from exc


    # Create new file for the final dataset
    final_data_fname = '../microSWIFT_data/final_dataset/' \
                       f'mission_{mission_num}.nc'

    # Open netCDF dataset
    rootgrp = nc.Dataset(final_data_fname, 'w',
                         clobber=True,
                         format='NETCDF4_CLASSIC')

    # Set overall file information
    rootgrp.title = f'DUNEX microSWIFT drifter - Mission {mission_num}'
    rootgrp.institution = 'University of Washington - Applied Physics Lab'
    rootgrp.source = 'Observations from microSWIFT drifters deployed ' \
                     'in the DUring Nearshore Events eXperiment (DUNEX)'
    rootgrp.Conventions = 'CF-1.6'
    rootgrp.Metadata_Conventions = 'Unidata Dataset Discovery v1.0'
    rootgrp.creator_country = 'USA'
    rootgrp.creator_email = 'erainvil@uw.edu'
    rootgrp.creator_name = 'EJ Rainville, Jim Thomson, Melissa Moulton, and ' \
                           'Morteza Derakhti at University of Washington - ' \
                           'Applied Physics Lab'
    rootgrp.creator_phone = '(303) 653-1226'
    rootgrp.creator_sector = 'academic'
    rootgrp.creator_state = 'Washington'
    rootgrp.featureType = 'trajectory'
    rootgrp.cdm_data_type = 'Trajectory'
    rootgrp.platform = 'microSWIFT wave buoy'
    rootgrp.publisher_country = 'USA'
    rootgrp.publisher_email = 'frfwebmaster@usace.army.mil'
    rootgrp.publisher_name = 'USACE/CHL/COAB'
    rootgrp.history = str(datetime.datetime.utcnow()) + ' Python'
    rootgrp.references = 'https://github.com/SASlabgroup/microSWIFT and ' \
                         'https://github.com/SASlabgroup/DUNEXMainExp'

    # Set Mission Specific Metadata
    mission_metadata(rootgrp, mission_num)

    # Create the dimensions for the dataset
    rootgrp.createDimension('time', cleaned_dataset['time'][:].shape[0])
    mission_time_nc = rootgrp.createVariable('time', np.float64, ('time',))
    mission_time_nc.units = cleaned_dataset['time'].units
    mission_time_nc.calendar = cleaned_dataset['time'].calendar
    mission_time_nc.standard_name = 'time'
    mission_time_nc.long_name = 'time'
    mission_time_nc[:] = cleaned_dataset['time'][:]

    microSWIFTs_on_mission = list(cleaned_dataset.groups.keys())
    microSWIFT_nums = [int(microSWIFT[11:]) for microSWIFT in microSWIFTs_on_mission]
    rootgrp.createDimension('trajectory', len(microSWIFTs_on_mission))
    trajectory_nc = rootgrp.createVariable('trajectory', np.intc, ('trajectory',))
    trajectory_nc.cf_role = 'trajectory_id'
    trajectory_nc.long_name = 'trajectory name'
    trajectory_nc.description = 'microSWIFT drift track ID. The ID number is ' \
                                'the same as the ID number of the microSWIFT ' \
                                'wave buoy.'
    trajectory_nc[:] = microSWIFT_nums

    # Add sampling frequency information to the netcdf file
    rootgrp.createDimension('single_value', 1)
    gps_freq = rootgrp.createVariable('gps_freq', np.intc, ('single_value'))
    gps_freq.long_name = 'GPS module sampling frequency'
    gps_freq.units = 's-1'
    gps_freq[:] = cleaned_dataset['gps_freq'][:]
    imu_freq = rootgrp.createVariable('imu_freq', np.intc, ('single_value'))
    imu_freq.long_name = 'IMU module sampling frequency'
    imu_freq.units = 's-1'
    imu_freq[:] = cleaned_dataset['imu_freq'][:]

    # Latitude
    lat_nc = rootgrp.createVariable('latitude',
                                    np.float64,
                                    ('trajectory', 'time'))
    lat_nc.standard_name = 'latitude'
    lat_nc.units = 'degrees_north'
    lat_nc.data_level = 'Level 1'
    lat_nc.coordinates = 'time latitude longitude'
    lat_nc.axis = 'Y'
    save_var_in_final_style('lat',
                            lat_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

    # Longitude
    lon_nc = rootgrp.createVariable('longitude',
                                    np.float64,
                                    ('trajectory', 'time'))
    lon_nc.standard_name = 'longitude'
    lon_nc.units = 'degrees_east'
    lon_nc.data_level = 'Level 1'
    lon_nc.coordinates = 'time latitude longitude'
    lon_nc.axis = 'X'
    save_var_in_final_style('lon',
                            lon_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

    # Body X Accelerations
    accel_x_body_nc = rootgrp.createVariable('acceleration_x_body',
                                             np.float64,
                                             ('trajectory', 'time'))
    accel_x_body_nc.long_name = 'X-axis acceleration in buoy reference frame'
    accel_x_body_nc.description = 'Acceleration across the major horizontal ' \
                                  'axis of the buoy in the reference frame ' \
                                  'of the buoy.'
    accel_x_body_nc.units = 'm s-2'
    accel_x_body_nc.data_level = 'Level 1'
    accel_x_body_nc.coordinates = 'time latitude longitude'
    save_var_in_final_style('accel_x_body',
                            accel_x_body_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

    # Body Y Accelerations
    accel_y_body_nc = rootgrp.createVariable('acceleration_y_body',
                                             np.float64,
                                             ('trajectory', 'time'))
    accel_y_body_nc.long_name = 'Y-axis acceleration in buoy reference frame'
    accel_y_body_nc.description = 'Acceleration across the minor horizontal ' \
                                  'axis of the buoy in the reference frame ' \
                                  'of the buoy.'
    accel_y_body_nc.units = 'm s-2'
    accel_y_body_nc.data_level = 'Level 1'
    accel_y_body_nc.coordinates = 'time latitude longitude'
    save_var_in_final_style('accel_y_body',
                            accel_y_body_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

    # Body Z Accelerations
    accel_z_body_nc = rootgrp.createVariable('acceleration_z_body',
                                             np.float64,
                                             ('trajectory', 'time'))
    accel_z_body_nc.long_name = 'Z-axis acceleration in buoy reference frame'
    accel_z_body_nc.description = 'Acceleration across the vertical ' \
                                  'axis of the buoy in the reference frame ' \
                                  'of the buoy (positive is down).'
    accel_z_body_nc.units = 'm s-2'
    accel_z_body_nc.data_level = 'Level 1'
    accel_z_body_nc.coordinates = 'time latitude longitude'
    save_var_in_final_style('accel_z_body',
                            accel_z_body_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

    # Gyro X Accelerations
    gyro_x_nc = rootgrp.createVariable('rotation_rate_x',
                                        np.float64,
                                        ('trajectory', 'time'))
    gyro_x_nc.long_name = 'X-axis gyroscope measurements in buoy ' \
                          'reference frame'
    gyro_x_nc.description = 'Rotation rate around the major horizontal ' \
                            'axis of the buoy.'
    gyro_x_nc.units = 'degrees s-1'
    gyro_x_nc.data_level = 'Level 1'
    gyro_x_nc.coordinates = 'time latitude longitude'
    save_var_in_final_style('gyro_x',
                            gyro_x_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

    # Gyro Y Accelerations
    gyro_y_nc = rootgrp.createVariable('rotation_rate_y',
                                        np.float64,
                                        ('trajectory', 'time'))
    gyro_y_nc.long_name = 'Y-axis gyroscope measurements in buoy ' \
                          'reference frame'
    gyro_y_nc.description = 'Rotation rate around the minor horizontal ' \
                            'axis of the buoy.'
    gyro_y_nc.units = 'degrees s-1'
    gyro_y_nc.data_level = 'Level 1'
    gyro_y_nc.coordinates = 'time latitude longitude'
    save_var_in_final_style('gyro_y',
                            gyro_y_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

    # Gyro Z Accelerations
    gyro_z_nc = rootgrp.createVariable('rotation_rate_z',
                                        np.float64,
                                        ('trajectory', 'time'))
    gyro_z_nc.long_name = 'Z-axis gyroscope measurements in buoy ' \
                          'reference frame'
    gyro_z_nc.description = 'Rotation rate around the vertical ' \
                            'axis of the buoy.'
    gyro_z_nc.units = 'degrees s-1'
    gyro_z_nc.data_level = 'Level 1'
    gyro_z_nc.coordinates = 'time latitude longitude'
    save_var_in_final_style('gyro_z',
                            gyro_z_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

    # Mag X Accelerations
    mag_x_nc = rootgrp.createVariable('magnetic_flux_density_x',
                                      np.float64,
                                      ('trajectory', 'time'))
    mag_x_nc.long_name = 'X-axis magnetic flux density measurements in buoy ' \
                          'reference frame'
    mag_x_nc.description = 'Magnetic flux density along the major horizontal ' \
                           'axis of the buoy.'
    mag_x_nc.units = 'uTesla'
    mag_x_nc.data_level = 'Level 1'
    mag_x_nc.coordinates = 'time latitude longitude'
    save_var_in_final_style('mag_x',
                            mag_x_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

    # Mag Y Accelerations
    mag_y_nc = rootgrp.createVariable('magnetic_flux_density_y',
                                      np.float64,
                                      ('trajectory', 'time'))
    mag_y_nc.long_name = 'Y-axis magnetic flux density measurements in buoy ' \
                          'reference frame'
    mag_y_nc.description = 'Magnetic flux density along the minor horizontal ' \
                           'axis of the buoy.'
    mag_y_nc.units = 'uTesla'
    mag_y_nc.data_level = 'Level 1'
    mag_y_nc.coordinates = 'time latitude longitude'
    save_var_in_final_style('mag_y',
                            mag_y_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

    # Mag Z Accelerations
    mag_z_nc = rootgrp.createVariable('magnetic_flux_density_z',
                                      np.float64,
                                      ('trajectory', 'time'))
    mag_z_nc.long_name = 'Z-axis magnetic flux density measurements in buoy ' \
                          'reference frame'
    mag_z_nc.description = 'Magnetic flux density along the major horizontal ' \
                           'axis of the buoy.'
    mag_z_nc.units = 'uTesla'
    mag_z_nc.data_level = 'Level 1'
    mag_z_nc.coordinates = 'time latitude longitude'
    save_var_in_final_style('mag_z',
                            mag_z_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

    # North-South Accelerations
    accel_ns_nc = rootgrp.createVariable('acceleration_ns',
                                         np.float64,
                                         ('trajectory', 'time'))
    accel_ns_nc.long_name = 'North-South acceleration in Earth ' \
                                 'reference frame'
    accel_ns_nc.description = 'Acceleration in the North-South direction in ' \
                              'the Earth reference frame. This is a ' \
                              'Level 2 data product that has been rotated ' \
                              'by the MATLAB AHRS indirect Kalman filter ' \
                              '(https://www.mathworks.com/help/fusion/ref/' \
                              'ahrsfilter-system-object.html) ' \
                              'and despiked using a PCHIP interpolation scheme' \
                              '(https://www.mathworks.com/help/matlab/ref/' \
                              'filloutliers.html).'
    accel_ns_nc.units = 'm s-2'
    accel_ns_nc.data_level = 'Level 2'
    accel_ns_nc.coordinates = 'time latitude longitude'
    save_var_in_final_style('accel_x',
                            accel_ns_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

    # East-West Accelerations
    accel_ew_nc = rootgrp.createVariable('acceleration_ew',
                                         np.float64,
                                         ('trajectory', 'time'))
    accel_ew_nc.long_name = 'East-West acceleration in Earth reference frame'
    accel_ew_nc.description = 'Acceleration in the East-West direction in ' \
                              'the Earth reference frame. This is a ' \
                              'Level 2 data product that has been rotated ' \
                              'by the MATLAB AHRS indirect Kalman filter ' \
                              '(https://www.mathworks.com/help/fusion/ref/' \
                              'ahrsfilter-system-object.html) ' \
                              'and despiked using a PCHIP interpolation scheme' \
                              '(https://www.mathworks.com/help/matlab/ref/' \
                              'filloutliers.html).'
    accel_ew_nc.units = 'm s-2'
    accel_ew_nc.data_level = 'Level 2'
    accel_ew_nc.coordinates = 'time latitude longitude'
    save_var_in_final_style('accel_y',
                            accel_ew_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

     # Up-Down Accelerations
    accel_ud_nc = rootgrp.createVariable('acceleration_ud',
                                         np.float64,
                                         ('trajectory', 'time'))
    accel_ud_nc.long_name = 'Up-Down acceleration in Earth reference frame'
    accel_ud_nc.description = 'Acceleration in the Up-Down direction in ' \
                              'the Earth reference frame. This is a ' \
                              'Level 2 data product that has been rotated ' \
                              'by the MATLAB AHRS indirect Kalman filter ' \
                              '(https://www.mathworks.com/help/fusion/ref/' \
                              'ahrsfilter-system-object.html) ' \
                              'and despiked using a PCHIP interpolation scheme' \
                              '(https://www.mathworks.com/help/matlab/ref/' \
                              'filloutliers.html). Down is Positive.'
    accel_ud_nc.units = 'm s-2'
    accel_ud_nc.data_level = 'Level 2'
    accel_ud_nc.coordinates = 'time latitude longitude'
    save_var_in_final_style('accel_z',
                            accel_ud_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)


    # East-West Velocity from GPS
    vel_ew_nc = rootgrp.createVariable('velocity_ew',
                                       np.float64,
                                       ('trajectory', 'time'))
    vel_ew_nc.long_name = 'East-West Velocity in Earth reference frame'
    vel_ew_nc.description = 'Velocity in the East-West direction ' \
                            'measured from the GPS module. This is a ' \
                            'Level 1 data product that has been despiked ' \
                            'and linearly interpolated onto the time ' \
                            'dimension.'
    vel_ew_nc.units = 'm s-1'
    vel_ew_nc.data_level = 'Level 1'
    vel_ew_nc.coordinates = 'time latitude longitude'
    save_var_in_final_style('u',
                            vel_ew_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

    # North-South Velocity from GPS
    vel_ns_nc = rootgrp.createVariable('velocity_ns',
                                       np.float64,
                                       ('trajectory', 'time'))
    vel_ns_nc.long_name = 'North-South Velocity in Earth reference frame'
    vel_ns_nc.description = 'Velocity in the North-South direction ' \
                            'measured from the GPS module. This is a ' \
                            'Level 1 data product that has been despiked ' \
                            'and linearly interpolated onto the time ' \
                            'dimension.'
    vel_ns_nc.units = 'm s-1'
    vel_ns_nc.data_level = 'Level 1'
    vel_ns_nc.coordinates = 'time latitude longitude'
    save_var_in_final_style('v',
                            vel_ns_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

    # Up-Down Velocity from IMU
    vel_ud_nc = rootgrp.createVariable('velocity_ud',
                                       np.float64,
                                       ('trajectory', 'time'))
    vel_ud_nc.long_name = 'Up-Down Velocity in Earth reference frame'
    vel_ud_nc.description = 'Velocity in the Up-Down direction integrated ' \
                            'from the Up-Down acceleration in the Earth ' \
                            'reference frame. This is a Level 2 data product ' \
                            'that has been integrated, filtered and despiked. '
    vel_ud_nc.units = 'm s-1'
    vel_ud_nc.data_level = 'Level 2'
    vel_ud_nc.coordinates = 'time latitude longitude'
    save_var_in_final_style('w',
                            vel_ud_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

    # X FRF location
    x_frf_nc = rootgrp.createVariable('xFRF',
                                      np.float64,
                                      ('trajectory', 'time'))
    x_frf_nc.long_name = 'x-coordinate in Local FRF Cartesian system'
    x_frf_nc.description = 'This is the cross-shore coordinate in the local ' \
                           'FRF coordinate system. Transformed from the ' \
                           'measured GPS locations.'
    x_frf_nc.units = 'm'
    x_frf_nc.data_level = 'Level 1'
    x_frf_nc.coordinates = 'time latitude longitude'
    save_var_in_final_style('xFRF',
                            x_frf_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

    # Y FRF location
    y_frf_nc = rootgrp.createVariable('yFRF',
                                      np.float64,
                                      ('trajectory', 'time'))
    y_frf_nc.long_name = 'y-coordinate in Local FRF Cartesian system'
    y_frf_nc.description = 'This is the along-shore coordinate in the local ' \
                           'FRF coordinate system. Transformed from the ' \
                           'measured GPS locations.'
    y_frf_nc.units = 'm'
    y_frf_nc.data_level = 'Level 1'
    y_frf_nc.coordinates = 'time latitude longitude'
    save_var_in_final_style('yFRF',
                            y_frf_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

    # Sea Surface Elevation
    sea_surf_elv_nc = rootgrp.createVariable('sea_surface_elevation',
                                             np.float64,
                                             ('trajectory', 'time'))
    sea_surf_elv_nc.long_name = 'Instantaneous sea surface elevation'
    sea_surf_elv_nc.description = 'This is the computed instantaneous sea ' \
                                  'surface elevation. This is a highly ' \
                                  'computed Level 2 data product. This is ' \
                                  'computed from the rotated Up-Down ' \
                                  'acceleration, double-integrated, filtered, ' \
                                  'and despiked.'
    sea_surf_elv_nc.units = 'm'
    sea_surf_elv_nc.data_level = 'Level 2'
    sea_surf_elv_nc.coordinates = 'time latitude longitude'
    save_var_in_final_style('eta',
                            sea_surf_elv_nc,
                            microSWIFTs_on_mission,
                            cleaned_dataset)

    rootgrp.close()


def mission_metadata(rootgrp, mission_num):
    """
    Set mission specific metadata from the DUNEX notes spreadsheet.
    """
    dunex_xlsx = pd.read_excel('../DUNEXMainExp_notes.xlsx', 'Main Experiment')
    deployer_1 = dunex_xlsx['Deployer 1'].iloc[mission_num]
    deployer_2 = dunex_xlsx['Deployer 2'].iloc[mission_num]
    retriever_1 = dunex_xlsx['Retriever 1'].iloc[mission_num]
    retriever_2 = dunex_xlsx['Retriever 2/ Notetaker'].iloc[mission_num]
    rootgrp.contributor_names = (f'{deployer_1}, {deployer_2}, ' \
                                f'{retriever_1}, {retriever_2}')
    rootgrp.contributor_role = 'Data Collectors'
    rootgrp.deployment_method = dunex_xlsx['Deployment Method'].iloc[mission_num]
    rootgrp.array_type = dunex_xlsx['Array Type'].iloc[mission_num]
    rootgrp.mission_start_time = dunex_xlsx['Start Time'].iloc[mission_num]
    rootgrp.mission_end_time = dunex_xlsx['End Time'].iloc[mission_num]
    rootgrp.deployment_notes = dunex_xlsx['Deployment Notes'].iloc[mission_num]

    # Get data cleaning notes
    # dunex_df = pd.read_excel('../DUNEXMainExp_notes.xlsx', 'microSWIFT Masks')
    # dunex_df['Mission Number'] = pd.Series(dunex_df['Mission Number']).fillna(method='ffill')
    # mission_masks = dunex_df[dunex_df['Mission Number'] == mission_num]
    # cleaning_notes = list(mission_masks[mission_masks['Mission Number'] == mission_num]['Notes'])
    # total_notes = ''
    # for note in cleaning_notes:
    #     if isinstance(note, str):
    #         total_notes += note
    #         total_notes += ' '
    #     else:
    #         continue
    # rootgrp.data_cleaning_notes = total_notes

def save_var_in_final_style(var_name,
                            nc_var,
                            microSWIFTs_on_mission,
                            cleaned_dataset):
    """
    Save the given variable in the new nc format from old file.

    Parameters
    ----------
    var_name : Name of the variable in the old netCDF
    """
    for n in range(len(microSWIFTs_on_mission)):
        nc_var[n,:] = cleaned_dataset[microSWIFTs_on_mission[n]][var_name][:]

if __name__ == '__main__':
    main()
