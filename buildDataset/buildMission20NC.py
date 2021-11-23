## buildMission20NC.py - Build Mission 20 Cleaned netCDF
# Import statements
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import netCDF4 as nc
from scipy import signal
from scipy import fft
from scipy import interpolate
from scipy import integrate
import cftime
import sys
import datetime

# Import DUNEX Tools
sys.path.append('..')
from tools import buildBasicMissionNC
from tools import missionView
from tools import microSWIFTTools

def main():
    '''
    @edwinrainville

    Description: Build and clean the mission 20 dataset.


    '''

    # Section 0: 
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

    # # Build a netCDF for the mission that has only raw time values read in 
    mission_nc_path = buildBasicMissionNC.main(mission_num=mission_num)

    # View the built mission netCDF
    missionView.main(mission_nc_path)

if __name__=='__main__':
    main()