# Build microSWIFT netCDF data structure from raw text data files
'''
@edwinrainville

'''
# Import statements
from netCDF4 import Dataset
import datetime
import glob
import numpy as np
import pynmea2
import cftime

# Define Project Directory 
project_dir = '/Volumes/DUNEXdata/DUNEXMainExp_Oct2021/'

# Define Data Directory
data_dir = 'microSWIFT_data/'

# Define Mission Number
# User input for mission number 
mission_num = int(input('Enter Mission Number: '))
mission_dir = 'mission_{}/'.format(mission_num)

# Define netCDF filename and path
ncfile_name = project_dir + data_dir + mission_dir + 'mission_{}.nc'.format(mission_num)

# Open netCDF dataset
rootgrp = Dataset(ncfile_name, 'w')

# Define microSWIFT num
microSWIFT_dir_list = glob.glob(project_dir + data_dir + mission_dir + 'microSWIFT_*')

# Loop through each directory to access data from mission
for microSWIFT_dir in microSWIFT_dir_list:

    # Define microSWIFT number from directory name
    microSWIFT_num = int(microSWIFT_dir[-2:])

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

    # Sort each list based on time before saving so that each data point is in chronological order
    imu_time = np.array(imu_time)
    imu_time_sorted_inds = imu_time.argsort()
    imu_time_sorted = imu_time[imu_time_sorted_inds]
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
    imu_time_dim = imugrp.createDimension('time', len(imu_time_sorted))

    # IMU Time Variable
    imu_time_nc = imugrp.createVariable('time', 'f8', ('time',))
    imu_time_nc.units = "hours since 0001-01-01 00:00:00.0"
    imu_time_nc.calendar = "gregorian"
    imu_time_num = cftime.date2num(imu_time_sorted, units=imu_time_nc.units,calendar=imu_time_nc.calendar)
    imu_time_nc[:] = imu_time_num

    # Accelerations
    accel_x_nc = imugrp.createVariable('accel_x', 'f8', ('time',))
    accel_x_nc[:] = accel_x_sorted
    accel_y_nc = imugrp.createVariable('accel_y', 'f8', ('time',))
    accel_y_nc[:] = accel_y_sorted
    accel_z_nc = imugrp.createVariable('accel_z', 'f8', ('time',))
    accel_z_nc[:] = accel_z_sorted

    # Magnetometer
    mag_x_nc = imugrp.createVariable('mag_x', 'f8', ('time',))
    mag_x_nc[:] = mag_x_sorted
    mag_y_nc = imugrp.createVariable('mag_y', 'f8', ('time',))
    mag_y_nc[:] = mag_y_sorted
    mag_z_nc = imugrp.createVariable('mag_z', 'f8', ('time',))
    mag_z_nc[:] = mag_z_sorted

    # Gyroscope
    gyro_x_nc = imugrp.createVariable('gyro_x', 'f8', ('time',))
    gyro_x_nc[:] = gyro_x_sorted
    gyro_y_nc = imugrp.createVariable('gyro_y', 'f8', ('time',))
    gyro_y_nc[:] = gyro_y_sorted
    gyro_z_nc = imugrp.createVariable('gyro_z', 'f8', ('time',))
    gyro_z_nc[:] = gyro_z_sorted

    # ------ GPS Data Read-in ------
    # Create GPS sub group
    gpsgrp = microSWIFTgroup.createGroup('GPS')

    # Get list of all GPS files
    gps_file_list = glob.glob(microSWIFT_dir + '/*GPS*')

    # Define lists for each variable
    gps_time = []
    lat = []
    lon = []
    u = []
    v = []
    z= []

    # Read in GPS data from each file in the GPS list
    for gps_file in gps_file_list:
        with open(gps_file, 'r') as file:
            
                for line in file:
                    if "GPGGA" in line:
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

                        # Read in other attributes from the GPGGA sentence
                        z.append(gpgga.altitude)
                        lat.append(gpgga.latitude)
                        lon.append(gpgga.longitude)
                    elif "GPVTG" in line:
                        if gpgga.gps_qual < 1:
                            continue
                        gpvtg = pynmea2.parse(line,check=True)   #grab gpvtg sentence
                        u.append(gpvtg.spd_over_grnd_kmph*np.cos(gpvtg.true_track)) #units are kmph
                        v.append(gpvtg.spd_over_grnd_kmph*np.sin(gpvtg.true_track)) #units are kmph
                    else: #if not GPGGA or GPVTG, continue to start of loop
                        continue

    # Sort each list based on time before saving so that each data point is in chronological order
    gps_time = np.array(gps_time)
    gps_time_sorted_inds = gps_time.argsort()
    gps_time_sorted = gps_time[gps_time_sorted_inds]
    lat_sorted = np.array(lat)[gps_time_sorted_inds]
    lon_sorted = np.array(lon)[gps_time_sorted_inds]
    z_sorted = np.array(z)[gps_time_sorted_inds]

    # Save GPS data to netCDF file
    # Create gps time dimension
    gps_time_dim = gpsgrp.createDimension('time', len(gps_time_sorted))

    # GPS Time Variable
    gps_time_nc = gpsgrp.createVariable('time', 'f8', ('time',))
    gps_time_nc.units = "hours since 0001-01-01 00:00:00.0"
    gps_time_nc.calendar = "gregorian"
    gps_time_num = cftime.date2num(gps_time_sorted, units=gps_time_nc.units,calendar=gps_time_nc.calendar)
    gps_time_nc[:] = gps_time_num

    # Locations
    lat_nc = gpsgrp.createVariable('lat', 'f8', ('time',))
    lat_nc[:] = lat_sorted
    lon_nc = gpsgrp.createVariable('lon', 'f8', ('time',))
    lon_nc[:] = lon_sorted
    z_nc = gpsgrp.createVariable('z', 'f8', ('time',))
    z_nc[:] = z_sorted

    # GPS Velocity variable
    gps_velocity = gpsgrp.createDimension('gps_velocity', len(u))
    u_nc = gpsgrp.createVariable('u', 'f8', ('gps_velocity',))
    u_nc[:] = u
    v_nc = gpsgrp.createVariable('v', 'f8', ('gps_velocity',))
    v_nc[:] = v







