## Building Mission 20 Data
import sys
sys.path.append('/Volumes/DUNEX/DUNEXMainExp/')

# Import Dataset tools
from tools.datasetTools import buildBasicMissionNC

def main():
    '''
    @edwinrainville

    Description: This script builds and cleans the dataset for microSWIFT mission 20 from the DUNEX Main
    experiment.

    '''

    # 0: Define basic values 
    mission_num = 20
    start_time = 'start'
    end_time = 'end'
   
    buildBasicMissionNC.mission_nc() 

    # 1: Read in Raw microSWIFT data and build initial netCDF file
    # mission_data_io(mission_num)

    # 2: Automated Cleaning of Dataset
    # microSWIFT_in_water(mission_num)
    # despike_microSWIFT_data(mission_num)

    # 3: Manual Cleaning of Dataset with highly detailed notes

    # 4: Compute Sea Surface Height from IMU and GPS data - This is an open question on how this will work 
    # compute_sea_surface_elevation(mission_num)

    # 5: Compute Wave Parameters - Hs, Tp, Dp, E(f) and E(f,theta) - This is an open question on how this will work 

    # 6: Find all wave breaking locations - This is based on Adam Browns Paper in 2018 - may need to be adjusted 

    # 7: Make spatial grid and organize all data into each grid cell spatially and in time - A grid will need to be decided

    # 8: Save all cleaned data into netCDF structure
    

# Run the Script
if __name__=='__main__':
    main()