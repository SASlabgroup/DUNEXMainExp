# Build microSWIFT netCDF data structure from raw text data files
# Import statements
import netCDF4 as nc
from scipy import interpolate
import datetime
import glob
import numpy as np
import pynmea2
import pandas as pd
import sys
import matlab.engine
import matlab

# Start the matlab engine
eng = matlab.engine.start_matlab()

# Import local Modules
sys.path.append('..')
from tools import microSWIFTTools

def main(mission_num=None):
    '''
    @edwinrainville

    TODO: Finish writing docstring
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
    dunex_xlsx = pd.read_excel(metadata_filename, 'Main Experiment')

    # microSWIFT Masks
    dunex_df = pd.read_excel(metadata_filename, 'microSWIFT Masks')
    dunex_df['Mission Number'] = pd.Series(dunex_df['Mission Number']).fillna(method='ffill')
    mission_masks = dunex_df[dunex_df['Mission Number'] == mission_num]

    # Define Mission Number
    # User input for mission number if not provided intially 
    if mission_num == 0:
        mission_num = int(input('Enter Mission Number: '))
    
    # Define Mission Directory 
    mission_dir = 'mission_{}/'.format(mission_num)

    # Define cleaned dataset directory
    cleaned_data_dir = 'cleanedDataset/'

    # Get start and end times
    start_time = datetime.datetime.fromisoformat(dunex_xlsx['Start Time'].iloc[mission_num])
    end_time = datetime.datetime.fromisoformat(dunex_xlsx['End Time'].iloc[mission_num])

    # Define netCDF filename and path
    ncfile_name = project_dir + data_dir + cleaned_data_dir + 'mission_{}.nc'.format(mission_num)

    # Open netCDF dataset
    rootgrp = nc.Dataset(ncfile_name, 'w', clobber=True)

    # Add sampling frequency information to the netcdf file
    single_val_dim = rootgrp.createDimension('single_value', 1)
    gps_freq = rootgrp.createVariable('gps_freq', 'f8', ('single_value'))
    gps_freq[:] = 4
    gps_freq.units = 'Hz'
    imu_freq = rootgrp.createVariable('imu_freq', 'f8', ('single_value'))
    imu_freq[:] = 12
    imu_freq.units = 'Hz'

    # Make mission time array and save to root group
    imu_time_step = datetime.timedelta(seconds=(1/imu_freq[0]))
    mission_time = np.arange(start_time, end_time, imu_time_step).astype(datetime.datetime)
    mission_time_dim = rootgrp.createDimension('time', len(mission_time))
    mission_time_nc = rootgrp.createVariable('time', 'f8', ('time',))
    mission_time_nc.units = "hours since 1970-01-01 00:00:00"
    mission_time_nc.calendar = "gregorian"
    mission_time_num = nc.date2num(mission_time, units=mission_time_nc.units,calendar=mission_time_nc.calendar)
    mission_time_nc[:] = mission_time_num

    # Define microSWIFT num
    microSWIFT_dir_list = glob.glob(project_dir + data_dir + mission_dir + 'microSWIFT_*')

    # Load in Bathymetry and Water Level Data to use to clean the dataset 
    bathy = nc.Dataset('../microSWIFT_data/FRFdata/FRF_geomorphology_DEMs_surveyDEM_20211021.nc')
    bathy_xFRF = bathy['xFRF'][:]
    bathy_yFRF = bathy['yFRF'][:]
    xFRF_grid, yFRF_grid = np.meshgrid(bathy_xFRF, bathy_yFRF)
    elevation = bathy['elevation'][0,:,:]

    # Load in water level data
    # url = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/oceanography/waterlevel/eopNoaaTide/2021/FRF-ocean_waterlevel_eopNoaaTide_202110.nc'
    waterLevel_dataset = nc.Dataset('../microSWIFT_data/FRFdata/FRF-ocean_waterlevel_eopNoaaTide_202110.nc')

    # Loop through each directory to access data from mission
    for microSWIFT_dir in microSWIFT_dir_list:
    
        # Define microSWIFT number from directory name
        microSWIFT_num = microSWIFT_dir[-2:]
        if microSWIFT_num[0] == '_':
            microSWIFT_num = int(microSWIFT_num[1:])
        else:
            microSWIFT_num = int(microSWIFT_num)

        # Check if there are both IMU and GPS files in the microSWIFT directory
        imu_file_list = glob.glob(microSWIFT_dir + '/*IMU*.dat')
        gps_file_list = glob.glob(microSWIFT_dir + '/*GPS*.dat')
        if (len(imu_file_list) > 0) and (len(gps_file_list) > 0):

            # ------ IMU Data Read-in ------
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
                with open(file, encoding="utf8", errors='ignore') as f:
                    lines = f.readlines()
                    # Line Structure: timestamp, accel_x, accel_y, accel_z, mag_x, mag_y, mag_z, gyro_x, gyro_y, gyro_z
                    for line in lines:
                        values = line.split(',')
                        if len(values) == 10:
                            # Read-In values from the IMU data files
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

            # Sort the time values 
            imu_time_sorted_inds = imu_time.argsort()

            # Sort the time and data based on sorting the time indices
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

            # Add in fraction of second of sampling frequency to each time value
            imu_time_with_millisecond = imu_time_sorted.copy()
            i = 1
            for n in np.arange(1, len(imu_time_sorted)):
                if imu_time_sorted[n] == imu_time_sorted[n-1]: 
                    imu_time_with_millisecond[n] = imu_time_sorted[n] + (i * imu_time_step) 
                    i += 1
                else:
                    # Restart i at one so that it can start again on the next second 
                    i = 1 

            # Sort the time to be within the mission time window from the notes spreadsheet
            current_ind = 0
            inds_in_mission = []
            for time in imu_time_with_millisecond:
                if time >= start_time and time <= end_time:
                    inds_in_mission.append(current_ind)
                    current_ind += 1
                else:
                    current_ind += 1

            # Check if there are any points within this time frame before trying to save them 
            if len(inds_in_mission) > 0:

                # Sort Indices to be within the mission time frame using the indices sorted above
                imu_time_sorted_in_mission = imu_time_with_millisecond[inds_in_mission]
                accel_x_sorted_in_mission = accel_x_sorted[inds_in_mission]
                accel_y_sorted_in_mission = accel_y_sorted[inds_in_mission]
                accel_z_sorted_in_mission = accel_z_sorted[inds_in_mission]
                mag_x_sorted_in_mission = mag_x_sorted[inds_in_mission]    
                mag_y_sorted_in_mission = mag_y_sorted[inds_in_mission]    
                mag_z_sorted_in_mission = mag_z_sorted[inds_in_mission]    
                gyro_x_sorted_in_mission = gyro_x_sorted[inds_in_mission]    
                gyro_y_sorted_in_mission = gyro_y_sorted[inds_in_mission]    
                gyro_z_sorted_in_mission = gyro_z_sorted[inds_in_mission]    

                # Map each index in the imu time series to an index in the overall time series
                mission_time_index_for_imu_value = np.searchsorted(mission_time, imu_time_sorted_in_mission, side='right')

                # If the last indexis the length of the array, replace it with the last index of the array
                mission_time_index_for_imu_value[mission_time_index_for_imu_value >= len(mission_time)] = -1
                
                # Make NaN vectors on mission_time for each variable and fill with values from the actual measurements
                # Accelerations
                accel_x_mission = np.nan * np.ones(len(mission_time))
                accel_x_mission[mission_time_index_for_imu_value] = accel_x_sorted_in_mission            
                accel_y_mission = np.nan * np.ones(len(mission_time))
                accel_y_mission[mission_time_index_for_imu_value] = accel_y_sorted_in_mission
                accel_z_mission = np.nan * np.ones(len(mission_time))
                accel_z_mission[mission_time_index_for_imu_value] = accel_z_sorted_in_mission

                # Magnetometer
                mag_x_mission = np.nan * np.ones(len(mission_time))
                mag_x_mission[mission_time_index_for_imu_value] = mag_x_sorted_in_mission
                mag_y_mission = np.nan * np.ones(len(mission_time))
                mag_y_mission[mission_time_index_for_imu_value] = mag_y_sorted_in_mission
                mag_z_mission = np.nan * np.ones(len(mission_time))
                mag_z_mission[mission_time_index_for_imu_value] = mag_z_sorted_in_mission

                # Gyroscope
                gyro_x_mission = np.nan * np.ones(len(mission_time))
                gyro_x_mission[mission_time_index_for_imu_value] = gyro_x_sorted_in_mission
                gyro_y_mission = np.nan * np.ones(len(mission_time))
                gyro_y_mission[mission_time_index_for_imu_value] = gyro_y_sorted_in_mission
                gyro_z_mission = np.nan * np.ones(len(mission_time))
                gyro_z_mission[mission_time_index_for_imu_value] = gyro_z_sorted_in_mission

                # Fill all gaps with linear interpolation 
                imu_nans = np.isnan(accel_x_mission)
                # accelerations
                accel_x_mission[imu_nans] = np.interp(mission_time_num[imu_nans], mission_time_num[~imu_nans], accel_x_mission[~imu_nans])
                accel_y_mission[imu_nans] = np.interp(mission_time_num[imu_nans], mission_time_num[~imu_nans], accel_y_mission[~imu_nans])
                accel_z_mission[imu_nans] = np.interp(mission_time_num[imu_nans], mission_time_num[~imu_nans], accel_z_mission[~imu_nans])
                # Gyroscope 
                gyro_x_mission[imu_nans] = np.interp(mission_time_num[imu_nans], mission_time_num[~imu_nans], gyro_x_mission[~imu_nans])
                gyro_y_mission[imu_nans] = np.interp(mission_time_num[imu_nans], mission_time_num[~imu_nans], gyro_y_mission[~imu_nans])
                gyro_z_mission[imu_nans] = np.interp(mission_time_num[imu_nans], mission_time_num[~imu_nans], gyro_z_mission[~imu_nans])
                # Magnetometer
                mag_x_mission[imu_nans] = np.interp(mission_time_num[imu_nans], mission_time_num[~imu_nans], mag_x_mission[~imu_nans])
                mag_y_mission[imu_nans] = np.interp(mission_time_num[imu_nans], mission_time_num[~imu_nans], mag_y_mission[~imu_nans])
                mag_z_mission[imu_nans] = np.interp(mission_time_num[imu_nans], mission_time_num[~imu_nans], mag_z_mission[~imu_nans])
                
                # Despike the acceleration, gyroscope and magnetometer measurements
                accel_x_despike = eng.filloutliers(matlab.double(accel_x_mission.tolist()), 'pchip', nargout=1)
                accel_y_despike = eng.filloutliers(matlab.double(accel_y_mission.tolist()), 'pchip', nargout=1)
                accel_z_despike = eng.filloutliers(matlab.double(accel_z_mission.tolist()), 'pchip', nargout=1)
                gyro_x_despike = eng.filloutliers(matlab.double(gyro_x_mission.tolist()), 'pchip', nargout=1)
                gyro_y_despike = eng.filloutliers(matlab.double(gyro_y_mission.tolist()), 'pchip', nargout=1)
                gyro_z_despike = eng.filloutliers(matlab.double(gyro_z_mission.tolist()), 'pchip', nargout=1)
                mag_x_despike = eng.filloutliers(matlab.double(mag_x_mission.tolist()), 'pchip', nargout=1)
                mag_y_despike = eng.filloutliers(matlab.double(mag_y_mission.tolist()), 'pchip', nargout=1)
                mag_z_despike = eng.filloutliers(matlab.double(mag_z_mission.tolist()), 'pchip', nargout=1)

                # Use MATLAB AHRS filter function to correct the accelerations to the earth frame of reference 
                accel_x_earth, accel_y_earth, accel_z_earth = eng.AHRSAccelCorrection(accel_x_despike, accel_y_despike, accel_z_despike, gyro_x_despike, gyro_y_despike, gyro_z_despike, mag_x_despike, mag_y_despike, mag_z_despike, nargout=3)

                # Convert all corrected accelerations back to numpy arrays
                accel_x_earth = np.squeeze(np.array(accel_x_earth))
                accel_y_earth = np.squeeze(np.array(accel_y_earth))
                accel_z_earth = np.squeeze(np.array(accel_z_earth))
                accel_x_despike = np.squeeze(np.array(accel_x_despike))
                accel_y_despike = np.squeeze(np.array(accel_y_despike))
                accel_z_despike = np.squeeze(np.array(accel_z_despike))

                # Integrate the vertical acceleration to get the vertical velocity and sea surface elevation
                vel_z_earth, z_earth = microSWIFTTools.computeEta(accel_z_earth)

                # Despike vertical velocity and sea surface elevation 
                vel_z_earth = np.squeeze(np.array(eng.filloutliers(matlab.double(vel_z_earth.tolist()), 'pchip', nargout=1)))
                z_earth = np.squeeze(np.array(eng.filloutliers(matlab.double(z_earth.tolist()), 'pchip', nargout=1)))

            # If there isn't any points within the mission - skip it
            else:
                continue

            # ------ GPS Data Read-in ------
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
                with open(gps_file, encoding="utf8", errors='ignore') as file:
                    for line in file:
                        if "GPGGA" in line:
                            linenum += 1
                            #check to see if we have lost GPS fix
                            try:
                                gpgga = pynmea2.parse(line)   #grab gpgga sentence and parse
                                if gpgga.gps_qual < 1:
                                    continue
                            except:
                                continue
                            else:
                            
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
                            try:
                                gpvtg = pynmea2.parse(line)   #grab gpvtg sentence
                                if type(gpvtg.spd_over_grnd_kmph) == float:
                                    vel_linenum.append(linenum)
                                    u.append(gpvtg.spd_over_grnd_kmph*np.cos(np.deg2rad(gpvtg.true_track)) * 0.277) # converting to m/s from km/hr
                                    v.append(gpvtg.spd_over_grnd_kmph*np.sin(np.deg2rad(gpvtg.true_track)) * 0.277) # converting to m/s from km/hr
                            except:
                                continue
                            
                        else: #if not GPGGA or GPVTG, continue to start of loop
                            linenum += 1
                            continue

            # Sort each list based on time before saving so that each data point is in chronological order
            gps_time = np.array(gps_time)
            gps_time_sorted_inds = gps_time.argsort()

            # Sorted GPS values
            gps_time_sorted = gps_time[gps_time_sorted_inds]
            gps_time_linenum = np.array(gps_time_linenum)[gps_time_sorted_inds]
            lat_sorted = np.array(lat)[gps_time_sorted_inds]
            lon_sorted = np.array(lon)[gps_time_sorted_inds]
            z_sorted = np.array(z)[gps_time_sorted_inds]

            # Sort the time to be within the mission time window from the notes spreadsheet
            current_ind = 0
            inds_in_mission = []
            for time in gps_time_sorted:
                if time >= start_time and time <= end_time:
                    inds_in_mission.append(current_ind)
                    current_ind += 1
                else:
                    current_ind += 1
                    continue
            
            if len(inds_in_mission) > 0:

                # Sort the values based on the indices within the mission time 
                gps_time_in_mission = np.array(gps_time_sorted)[inds_in_mission]
                gps_time_linenum_in_mission = gps_time_linenum[inds_in_mission]
                lat_sorted_in_mission = lat_sorted[inds_in_mission]
                lon_sorted_in_mission = lon_sorted[inds_in_mission]
                z_sorted_in_mission = z_sorted[inds_in_mission]

                # Interpolate each value onto the overall mission time
                gps_time_num = nc.date2num(gps_time_in_mission, units=mission_time_nc.units,calendar=mission_time_nc.calendar) 
                
                # Latitude interpolation and saving to netCDF in the microSWIFT group
                lat_interp_func = interpolate.interp1d(gps_time_num, lat_sorted_in_mission, bounds_error=False, fill_value='NaN')
                lat_interpolated = lat_interp_func(mission_time_num)

                # Longitude interpolation and saving to netCDF in the microSWIFT group
                lon_interp_func = interpolate.interp1d(gps_time_num, lon_sorted_in_mission, bounds_error=False, fill_value='NaN')
                lon_interpolated = lon_interp_func(mission_time_num)

                # GPS Elevation interpolation and saving to netCDF in the microSWIFT group
                z_interp_func = interpolate.interp1d(gps_time_num, z_sorted_in_mission, bounds_error=False, fill_value='NaN')
                z_interpolated = z_interp_func(mission_time_num)

                # Compute FRF x and y locations 
                x, y = microSWIFTTools.transform2FRF(lat=lat_interpolated, lon=lon_interpolated)
                
                # Interpolate times from GPGGA time to the GPVTG time
                extrap_func = interpolate.interp1d(gps_time_linenum_in_mission, gps_time_num, fill_value='extrapolate')
                gps_vel_time = extrap_func(vel_linenum)

                # Sort GPS velocites based on gps time
                gps_vel_time = np.array(gps_vel_time)
                gps_vel_time_sorted_inds = gps_vel_time.argsort()
                gps_vel_time_sorted = gps_vel_time[gps_vel_time_sorted_inds]
                u_sorted = np.array(u)[gps_vel_time_sorted_inds]
                v_sorted = np.array(v)[gps_vel_time_sorted_inds]

                # Interpolate the GPS velocities onto the mission time array
                gps_u_interp_func = interpolate.interp1d(gps_vel_time_sorted, u_sorted, bounds_error=False, fill_value='NaN')
                gps_u_mission = gps_u_interp_func(mission_time_num)
                gps_v_interp_func = interpolate.interp1d(gps_vel_time_sorted, v_sorted, bounds_error=False, fill_value='NaN')
                gps_v_mission = gps_v_interp_func(mission_time_num)
            
            # If there aren't any points in the mission time skip 
            else:
                continue
                
            # ------  Clean the dataset using the indices in the mission notesheet and bathymetry --------
            # Find water level closest to during mission 
            waterLevel_in_mission = np.interp(np.median(mission_time_num), waterLevel_dataset['time'][:], waterLevel_dataset['waterLevel'][:])
            
            # Adjust the bathymetry with the water level 
            depth = elevation + waterLevel_in_mission
            depth[depth < 0] = 0
            depth[depth > 0 ] = 1

            # Add a 2 meter buffer to the furtherest cross shore point
            buffer = 2 # units are meters

            # Find the max x where there is beach showing and add buffer
            max_x_row = []
            for n in np.arange(len(bathy_yFRF)):
                max_x_row.append(bathy_xFRF[np.where(depth[n,:] == 0)[0][0]])
            max_x = np.max(max_x_row)
            max_x += buffer

            if mission_num == 51 or mission_num == 52:
                max_x = -100 # for mission 51 and 52 since they have been cleaned very throoughyl with hand written notes

            # Find first index in the xFRF data that is less than the beach location
            beach_mask_ind = np.where(x <= max_x)[0].tolist()

            # Mask all indices where the microSWIFT is on the beach
            microSWIFT_variables = [accel_x_despike, accel_y_despike, accel_z_despike,
                                    accel_x_earth, accel_y_earth, accel_z_earth, 
                                    gps_u_mission, gps_v_mission, vel_z_earth, 
                                    x, y, z_earth, mag_x_mission, mag_y_mission, 
                                    mag_z_mission, gyro_x_mission, gyro_y_mission,
                                    gyro_z_mission, lat_interpolated, lon_interpolated]
            for variable in microSWIFT_variables:
                for n in beach_mask_ind:
                    variable[n] = np.NaN
            
            # Mask Individual indices from the Mission Note Sheet
            microSWIFTs_masked_on_mission = list(mission_masks[mission_masks['Mission Number'] == mission_num]['microSWIFT ID'])

            if 'microSWIFT_{}'.format(microSWIFT_num) in microSWIFTs_masked_on_mission:
                # Get the mask values from spreadsheet 
                individual_microSWIFT_mask = mission_masks[mission_masks['microSWIFT ID'] == 'microSWIFT_{}'.format(microSWIFT_num)]
                
                # Mask from Begininng to start mask end index 
                start_mask_end_index = int(individual_microSWIFT_mask['Start Mask End Index'].item())
                end_mask_start_index = int(individual_microSWIFT_mask['End Mask Start Index'].item())
                slice_list = []
                additional_masked_points = []
                if individual_microSWIFT_mask['Additional Masking indices'].item() != 0:
                    # sort between ranges we want to mask and individual points
                    for val in individual_microSWIFT_mask['Additional Masking indices'].item().split(','):
                        if ':' in val:
                            slice_list.append(slice(int(val.split(':')[0]), int(val.split(':')[1])))
                        else:
                            additional_masked_points.append(int(val))

                    # Convert additional points list to numpy array
                    additional_masked_points = np.array(additional_masked_points)

                # Mask all indices in the mask list
                for variable in microSWIFT_variables:
                    variable[:start_mask_end_index] = np.NaN
                    variable[end_mask_start_index:] = np.NaN
                    if len(additional_masked_points) > 0:
                        variable[additional_masked_points] = np.NaN
                    if len(slice_list) > 0:
                        for inds in slice_list:
                            variable[inds] = np.NaN

            # Save all loaded in data to a microSWIFT subgroup 
            # Check that this microSWIFT has all data before making a subgroup and saving
            if np.all(np.isnan(lat))==False and np.all(np.isnan(lon))==False and np.all(np.isnan(accel_x_despike))==False:
                if (np.count_nonzero(~np.isnan(accel_z_despike)) > 2000 and (np.nanmean(accel_z_despike) > 9) and (np.nanmean(accel_z_despike) < 11)) == True:
                    # if mission_num == 51 or mission_num == 52: 
                    # Create netcdf group for microSWIFT
                    microSWIFTgroup = rootgrp.createGroup('microSWIFT_{}'.format(microSWIFT_num))

                    # Save IMU data to microSWIFT group
                    # Acceleration - Body Frame of Reference
                    accel_x_nc = microSWIFTgroup.createVariable('accel_x_body', 'f8', ('time',))
                    accel_x_nc.units = 'm/s^2'
                    accel_x_nc[:] = accel_x_despike
                    accel_y_nc = microSWIFTgroup.createVariable('accel_y_body', 'f8', ('time',))
                    accel_y_nc.units = 'm/s^2'
                    accel_y_nc[:] = accel_y_despike
                    accel_z_nc = microSWIFTgroup.createVariable('accel_z_body', 'f8', ('time',))
                    accel_z_nc.units = 'm/s^2'
                    accel_z_nc[:] = accel_z_despike

                    # Acceleration - Earth Frame of Reference 
                    accel_x_nc_earth = microSWIFTgroup.createVariable('accel_x', 'f8', ('time',))
                    accel_x_nc_earth.units = 'm/s^2'
                    accel_x_nc_earth[:] = accel_x_earth
                    accel_y_nc_earth = microSWIFTgroup.createVariable('accel_y', 'f8', ('time',))
                    accel_y_nc_earth.units = 'm/s^2'
                    accel_y_nc_earth[:] = accel_y_earth
                    accel_z_nc_earth = microSWIFTgroup.createVariable('accel_z', 'f8', ('time',))
                    accel_z_nc_earth.units = 'm/s^2'
                    accel_z_nc_earth[:] = accel_z_earth

                    # Velocity
                    u_nc = microSWIFTgroup.createVariable('u', 'f8', ('time',))
                    u_nc.units = 'm/s'
                    u_nc[:] = gps_u_mission
                    v_nc = microSWIFTgroup.createVariable('v', 'f8', ('time',))
                    v_nc.units = 'm/s'
                    v_nc[:] = gps_v_mission
                    w_nc = microSWIFTgroup.createVariable('w', 'f8', ('time',))
                    w_nc.units = 'm/s'
                    w_nc[:] = vel_z_earth

                    # Positions
                    x_frf_nc = microSWIFTgroup.createVariable('xFRF', 'f8', ('time',))
                    x_frf_nc.units = 'meters'
                    x_frf_nc[:] = x
                    y_frf_nc = microSWIFTgroup.createVariable('yFRF', 'f8', ('time',))
                    y_frf_nc.units = 'meters'
                    y_frf_nc[:] = y
                    z_nc = microSWIFTgroup.createVariable('eta', 'f8', ('time',))
                    z_nc.units = 'meters'
                    z_nc[:] = z_earth

                    # Magnetometer
                    mag_x_nc = microSWIFTgroup.createVariable('mag_x', 'f8', ('time',))
                    mag_x_nc.units = 'uTeslas'
                    mag_x_nc[:] = mag_x_mission
                    mag_y_nc = microSWIFTgroup.createVariable('mag_y', 'f8', ('time',))
                    mag_y_nc.units = 'uTeslas'
                    mag_y_nc[:] = mag_y_mission
                    mag_z_nc = microSWIFTgroup.createVariable('mag_z', 'f8', ('time',))
                    mag_z_nc.units = 'uTeslas'
                    mag_z_nc[:] = mag_z_mission

                    # Gyroscope
                    gyro_x_nc = microSWIFTgroup.createVariable('gyro_x', 'f8', ('time',))
                    gyro_x_nc.units = 'degrees/sec'
                    gyro_x_nc[:] = gyro_x_mission
                    gyro_y_nc = microSWIFTgroup.createVariable('gyro_y', 'f8', ('time',))
                    gyro_y_nc.units = 'degrees/sec'
                    gyro_y_nc[:] = gyro_y_mission
                    gyro_z_nc = microSWIFTgroup.createVariable('gyro_z', 'f8', ('time',))
                    gyro_z_nc.units = 'degrees/sec'
                    gyro_z_nc[:] = gyro_z_mission

                    # Lat and Lon
                    lat_nc = microSWIFTgroup.createVariable('lat', 'f8', ('time',))
                    lat_nc.units = 'degrees_north'
                    lat_nc[:] = lat_interpolated
                    lon_nc = microSWIFTgroup.createVariable('lon', 'f8', ('time',))
                    lon_nc.units = 'degrees_east'
                    lon_nc[:] = lon_interpolated

    # Close the dataset
    rootgrp.close()

    # Return the name of the netCDF that was created
    return ncfile_name

    # Close the matlab engine
    eng.quit()

# Run the Script
if __name__=='__main__':
    main()
