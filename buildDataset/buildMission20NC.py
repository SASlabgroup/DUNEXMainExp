## Building Mission 20 Data
# Import Statements
import numpy as np
import pandas as pd
import datetime
import sys
sys.path.append('..')

# Import Dataset tools
from tools import buildBasicMissionNC
from tools import microSWIFTTools

def main():
    '''
    @edwinrainville

    Description: This script builds and cleans the dataset for microSWIFT mission 20 from the DUNEX Main
    experiment.

    '''
    # ---------- Step 0: Define basic values ------------
    mission_num = 20
    
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

    # Get start and end times
    start_time = datetime.datetime.fromisoformat(dunex_xlsx['Start Time'].iloc[mission_num])
    end_time = datetime.datetime.fromisoformat(dunex_xlsx['End Time'].iloc[mission_num])
   
    # ---------- Step 1: Build Basic netCDF file from mission microSWIFT data ------------
    # Build a netCDF for the mission that has only raw time values read in 
    buildBasicMissionNC.main(mission_num=mission_num)

    # Match GPS time steps 


    # Match IMU time steps

    # ---------- Step 2: Automated Cleaning of Dataset ------------
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