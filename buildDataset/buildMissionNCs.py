## buildMission20NC.py - Build Mission 20 Cleaned netCDF
# Import statements
import numpy as np
import pandas as pd
import sys
import datetime

# Import DUNEX Tools
sys.path.append('..')
from tools import buildBasicMissionNC

def main():
    '''
    @edwinrainville

    Description: Build and clean the basic mission datasets.
    TODO: Finish writing description
    '''
    # Define number of missions
    number_of_missions = 81
    
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

    # Loop through each mission and build netCDF files
    for mission_num in np.arange(20, number_of_missions):
        # Get start and end times
        start_time = datetime.datetime.fromisoformat(dunex_xlsx['Start Time'].iloc[mission_num])
        end_time = datetime.datetime.fromisoformat(dunex_xlsx['End Time'].iloc[mission_num])

        # Build a netCDF for the mission that has only raw time values read in 
        print('building mission {}'.format(mission_num))
        mission_nc_path = buildBasicMissionNC.main(mission_num=mission_num)

if __name__=='__main__':
    main()