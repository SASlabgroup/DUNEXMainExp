# Build microSWIFT netCDF data structure from raw text data files
# Import statements
from netCDF4 import Dataset
import netCDF4 as nc
from scipy import interpolate
import datetime
import glob
import numpy as np
import pynmea2
import cftime
import microSWIFTTools
import pandas as pd

def main():
    '''
    @edwinrainville

    Description: This script takes all the raw offloaded microSWIFT data and loads it in to a well organized and formatted 
                netCDF file that can then be easily analyzed for each mission.

    '''

    # Define Project Directory 
    project_dir = '../'

    # Define Data Directory
    data_dir = 'microSWIFT_data/'

    # Define Metadata Excel sheet name
    metadata_name = 'DUNEXMainExp_notes.xlsx'

    # Combine file name and project Directory
    metadata_filename = project_dir + metadata_name

    # Create dataframe object from DUNEX MetaData SpreadSheet
    dunex_xlsx = pd.read_excel(metadata_filename)

    # Define Mission Number
    # User input for mission number 
    mission_num = int(input('Enter Mission Number: '))
    mission_dir = 'mission_{}/'.format(mission_num)

    # Get start and end times
    start_time = datetime.datetime.fromisoformat(dunex_xlsx['Start Time'].iloc[mission_num])
    end_time = datetime.datetime.fromisoformat(dunex_xlsx['End Time'].iloc[mission_num])

    # Define netCDF filename and path
    ncfile_name = project_dir + data_dir + mission_dir + 'mission_{}.nc'.format(mission_num)

    # Open netCDF dataset
    rootgrp = Dataset(ncfile_name, 'w', clobber=True)

    # Add sampling frequency information to the netcdf file
    single_val_dim = rootgrp.createDimension('single_value', 1)
    gps_freq = rootgrp.createVariable('gps_freq', 'f8', ('single_value'))
    gps_freq[:] = 4
    gps_freq.units = 'Hz'
    imu_freq = rootgrp.createVariable('imu_freq', 'f8', ('single_value'))
    imu_freq[:] = 12
    imu_freq.units = 'Hz'

    # Define microSWIFT num
    microSWIFT_dir_list = glob.glob(project_dir + data_dir + mission_dir + 'microSWIFT_*')

    # Loop through each directory to access data from mission
    for microSWIFT_dir in microSWIFT_dir_list:

        # Define microSWIFT number from directory name
        microSWIFT_num = microSWIFT_dir[-2:]
        if microSWIFT_num[0] == '_':
            microSWIFT_num = int(microSWIFT_num[1:])
        else:
            microSWIFT_num = int(microSWIFT_num)
        

        print(microSWIFT_num)

        # Create netcdf group for microSWIFT
        microSWIFTgroup = rootgrp.createGroup('microSWIFT_{}'.format(microSWIFT_num))

        # ------ IMU Data Read-in ------
        # Create IMU sub group
        imugrp = microSWIFTgroup.createGroup('IMU')

        # Get list of all IMU files
        imu_file_list = glob.glob(microSWIFT_dir + '/*IMU*')

        # Define lists for each variable
        imu_time = []
        accel_x = []
        accel_y = []
        accel_z = []
        mag_x =[]
        mag_y = []
        mag_z = []
        gyro_x = []
        gyro_y = []
        gyro_z = []

        # Loop through each file and read in data from each line 
        for file in imu_file_list:
            print(file)
            with open(file) as f:
                lines = f.readlines()
                # Line Structure: timestamp, accel_x, accel_y, accel_z, mag_x, mag_y, mag_z, gyro_x, gyro_y, gyro_z
                for line in lines:
                    values = line.split(',')
                    if len(values) == 10:
                        imu_time.append(datetime.datetime.fromisoformat(values[0]))
                        accel_x.append(float(values[1]))
                        accel_y.append(float(values[2]))
                        accel_z.append(float(values[3]))
                        mag_x.append(float(values[4]))
                        mag_y.append(float(values[5]))
                        mag_z.append(float(values[6]))
                        gyro_x.append(float(values[7]))
                        gyro_y.append(float(values[8]))
                        gyro_z.append(float(values[9]))
                    else:
                        continue

        # Sort the time to be within the mission time window from the notes spreadsheet
        imu_time_in_mission = []
        for time in imu_time:
            if time >= start_time and time <= end_time:
                imu_time_in_mission.append(time)
            else:
                continue

        # Sort each list based on time before saving so that each data point is in chronological order
        imu_time_in_mission = np.array(imu_time_in_mission)
        imu_time_sorted_inds = imu_time_in_mission.argsort()
        imu_time_sorted = imu_time_in_mission[imu_time_sorted_inds]
        accel_x_sorted = np.array(accel_x)[imu_time_sorted_inds]
        accel_y_sorted = np.array(accel_y)[imu_time_sorted_inds]
        accel_z_sorted = np.array(accel_z)[imu_time_sorted_inds]
        mag_x_sorted = np.array(mag_x)[imu_time_sorted_inds]    
        mag_y_sorted = np.array(mag_y)[imu_time_sorted_inds]    
        mag_z_sorted = np.array(mag_z)[imu_time_sorted_inds]    
        gyro_x_sorted = np.array(gyro_x)[imu_time_sorted_inds]    
        gyro_y_sorted = np.array(gyro_y)[imu_time_sorted_inds]    
        gyro_z_sorted = np.array(gyro_z)[imu_time_sorted_inds]    

        # Create IMU dimensions and write data to netCDF file
        # Create imu time dimension
        imu_time_dim = imugrp.createDimension('time', imu_time_sorted.shape[0])

        # IMU Time Variable
        imu_time_nc = imugrp.createVariable('time', 'f8', ('time',))
        imu_time_nc.units = "hours since 1970-01-01 00:00:00"
        imu_time_nc.calendar = "gregorian"
        imu_time_num = nc.date2num(imu_time_sorted, units=imu_time_nc.units,calendar=imu_time_nc.calendar)
        imu_time_nc[:] = imu_time_num

        # Accelerations
        accel_x_nc = imugrp.createVariable('accel_x', 'f8', ('time',))
        accel_x_nc.units = 'm/s^2'
        accel_x_nc[:] = accel_x_sorted
        accel_y_nc = imugrp.createVariable('accel_y', 'f8', ('time',))
        accel_y_nc.units = 'm/s^2'
        accel_y_nc[:] = accel_y_sorted
        accel_z_nc = imugrp.createVariable('accel_z', 'f8', ('time',))
        accel_z_nc.units = 'm/s^2'
        accel_z_nc[:] = accel_z_sorted

        # Magnetometer
        mag_x_nc = imugrp.createVariable('mag_x', 'f8', ('time',))
        mag_x_nc.units = 'uTeslas'
        mag_x_nc[:] = mag_x_sorted
        mag_y_nc = imugrp.createVariable('mag_y', 'f8', ('time',))
        mag_y_nc.units = 'uTeslas'
        mag_y_nc[:] = mag_y_sorted
        mag_z_nc = imugrp.createVariable('mag_z', 'f8', ('time',))
        mag_z_nc.units = 'uTeslas'
        mag_z_nc[:] = mag_z_sorted

        # Gyroscope
        gyro_x_nc = imugrp.createVariable('gyro_x', 'f8', ('time',))
        gyro_x_nc.units = 'degrees/sec'
        gyro_x_nc[:] = gyro_x_sorted
        gyro_y_nc = imugrp.createVariable('gyro_y', 'f8', ('time',))
        gyro_y_nc.units = 'degrees/sec'
        gyro_y_nc[:] = gyro_y_sorted
        gyro_z_nc = imugrp.createVariable('gyro_z', 'f8', ('time',))
        gyro_z_nc.units = 'degrees/sec'
        gyro_z_nc[:] = gyro_z_sorted

        # ------ GPS Data Read-in ------
        # Create GPS sub group
        gpsgrp = microSWIFTgroup.createGroup('GPS')

        # Get list of all GPS files
        gps_file_list = glob.glob(microSWIFT_dir + '/*GPS*')

        # Define lists for each variable
        gps_time = []
        gps_time_linenum = []
        lat = []
        lon = []
        u = []
        v = []
        z= []
        vel_linenum = []
        linenum = 0

        # Read in GPS data from each file in the GPS list
        for gps_file in gps_file_list:
            print(gps_file)
            
            with open(gps_file, 'r') as file:
                
                    for line in file:
                        if "GPGGA" in line:
                            linenum += 1
                            gpgga = pynmea2.parse(line,check=True)   #grab gpgga sentence and parse
                            #check to see if we have lost GPS fix
                            if gpgga.gps_qual < 1:
                                continue
                            
                            # Create datetime from timestamp and file name
                            # Convert Month string to month num
                            month_str = gps_file[-21:-18]
                        
                            # January
                            if month_str == 'Jan':
                                month_num = '01'
                            # February
                            if month_str == 'Feb':
                                month_num = '02'
                            # March
                            if month_str == 'Mar':
                                month_num = '03'
                            # April
                            if month_str == 'Apr':
                                month_num = '04'
                            # May
                            if month_str == 'May':
                                month_num = '05'
                            # June
                            if month_str == 'Jun':
                                month_num = '06'
                            # July
                            if month_str == 'Jul':
                                month_num = '07'
                            # August
                            if month_str == 'Aug':
                                month_num = '08'
                            # September
                            if month_str == 'Sep':
                                month_num = '09'
                            # October 
                            if month_str == 'Oct':
                                month_num = '10'
                            # November
                            if month_str == 'Nov':
                                month_num = '11'
                            # December
                            if month_str == 'Dec':
                                month_num = '12'

                            # Compute Datetime
                            date_str = '{0}-{1}-{2}'.format(gps_file[-18:-14], month_num, gps_file[-23:-21])
                            gps_date = datetime.date.fromisoformat(date_str)
                            gps_datetime = datetime.datetime.combine(gps_date, gpgga.timestamp)
                            gps_time.append(gps_datetime)
                            gps_time_linenum_val = linenum
                            gps_time_linenum.append(gps_time_linenum_val)
                            # Read in other attributes from the GPGGA sentence
                            z.append(gpgga.altitude)
                            lat.append(gpgga.latitude)
                            lon.append(gpgga.longitude)
                        elif "GPVTG" in line:
                            linenum += 1
                            if gpgga.gps_qual < 1:
                                continue
                            vel_linenum.append(linenum)
                            gpvtg = pynmea2.parse(line,check=True)   #grab gpvtg sentence
                            u.append(gpvtg.spd_over_grnd_kmph*np.cos(gpvtg.true_track)) #units are kmph
                            v.append(gpvtg.spd_over_grnd_kmph*np.sin(gpvtg.true_track)) #units are kmph
                        else: #if not GPGGA or GPVTG, continue to start of loop
                            linenum += 1
                            continue

        # Sort the time to be within the mission time window from the notes spreadsheet
        gps_time_in_mission = []
        for time in gps_time:
            if time >= start_time and time <= end_time:
                gps_time_in_mission.append(time)
            else:
                continue

        # Sort each list based on time before saving so that each data point is in chronological order
        gps_time_in_mission = np.array(gps_time_in_mission)
        gps_time_sorted_inds = gps_time_in_mission.argsort()

        # Sorted GPS values
        gps_time_sorted = gps_time_in_mission[gps_time_sorted_inds]
        gps_time_linenum = np.array(gps_time_linenum)[gps_time_sorted_inds]
        lat_sorted = np.array(lat)[gps_time_sorted_inds]
        lon_sorted = np.array(lon)[gps_time_sorted_inds]
        z_sorted = np.array(z)[gps_time_sorted_inds]

        # Save GPS data to netCDF file
        gps_time_dim = gpsgrp.createDimension('time', len(gps_time_sorted))

        # GPS Time Variable
        gps_time_nc = gpsgrp.createVariable('time', 'f8', ('time',))
        gps_time_nc.units = "hours since 1970-01-01 00:00:00"
        gps_time_nc.calendar = "standard"
        gps_time_num = nc.date2num(gps_time_sorted, units=gps_time_nc.units,calendar=gps_time_nc.calendar)
        gps_time_nc[:] = gps_time_num

        # Locations
        lat_nc = gpsgrp.createVariable('lat', 'f8', ('time',))
        lat_nc.units = 'degrees_north'
        lat_nc[:] = lat_sorted
        lon_nc = gpsgrp.createVariable('lon', 'f8', ('time',))
        lon_nc.units = 'degrees_east'
        lon_nc[:] = lon_sorted
        z_nc = gpsgrp.createVariable('z', 'f8', ('time',))
        z_nc[:] = z_sorted

        # Compute FRF x and y locations 
        x, y = microSWIFTTools.transform2FRF(lat=lat_sorted, lon=lon_sorted)
        x_frf_nc = gpsgrp.createVariable('x_frf', 'f8', ('time',))
        x_frf_nc.units = 'meters'
        x_frf_nc[:] = x
        y_frf_nc = gpsgrp.createVariable('y_frf', 'f8', ('time',))
        y_frf_nc.units = 'meters'
        y_frf_nc[:] = y
        
        # Interpolate times from GPGGA time to the GPVTG time
        # gps_time_num_notsorted = cftime.date2num(gps_time, units=gps_time_nc.units,calendar=gps_time_nc.calendar)
        extrap_func = interpolate.interp1d(gps_time_linenum, gps_time_num, fill_value='extrapolate')
        gps_vel_time = extrap_func(vel_linenum)

        # Sort GPS velocites based on gps time
        gps_vel_time = np.array(gps_vel_time)
        gps_vel_time_sorted_inds = gps_vel_time.argsort()
        gps_vel_time_sorted = gps_vel_time[gps_vel_time_sorted_inds]
        u_sorted = np.array(u)[gps_vel_time_sorted_inds]
        v_sorted = np.array(v)[gps_vel_time_sorted_inds]

        # Add all GPS Velocity values to the mission netcdf file
        gps_vel_time_dim = gpsgrp.createDimension('gps_velocity_time', len(vel_linenum))
        gps_vel_time_nc = gpsgrp.createVariable('gps_vel_time', 'f8', ('gps_velocity_time',))
        gps_vel_time_nc.units = "hours since 1970-01-01 00:00:00"
        gps_vel_time_nc.calendar = "standard"
        gps_vel_time_nc[:] = gps_vel_time
        u_nc = gpsgrp.createVariable('u', 'f8', ('gps_velocity_time',))
        u_nc.units = 'm/s'
        u_nc[:] = u_sorted
        v_nc = gpsgrp.createVariable('v', 'f8', ('gps_velocity_time',))
        v_nc.units = 'm/s'
        v_nc[:] = v_sorted

# Run the Script
if __name__=='__main__':
    main()
