# Wave Conditions at the FRF 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import netCDF4 as nc
import matplotlib.patches as patches
import cftime

def main():
    # Create dataframe object from DUNEX MetaData SpreadSheet
    dunex_xlsx = pd.read_excel('../DUNEXMainExp_notes.xlsx')

    # Get start and end times from notes sheet
    start_times = []
    end_times = []
    for n in np.arange(len(dunex_xlsx['Start Time'])):
        start_times.append(datetime.datetime.fromisoformat(dunex_xlsx['Start Time'][n]))
        end_times.append(datetime.datetime.fromisoformat(dunex_xlsx['End Time'][n]))

    # Get AWAC data
    awac4p5m = nc.Dataset('../microSWIFT_data/FRFdata/FRF-ocean_waves_awac-4.5m_202110.nc')
    awac4p5m_time = cftime.num2pydate(awac4p5m['time'], units=awac4p5m['time'].units, calendar=awac4p5m['time'].calendar)
    awac6m = nc.Dataset('../microSWIFT_data/FRFdata/FRF-ocean_waves_awac-6m_202110.nc')
    awac6m_time = cftime.num2pydate(awac6m['time'], units=awac6m['time'].units, calendar=awac6m['time'].calendar)
    array8m = nc.Dataset('../microSWIFT_data/FRFdata/FRF-ocean_waves_8m-array_202110.nc')
    array8m_time = cftime.num2pydate(array8m['time'], units=array8m['time'].units, calendar=array8m['time'].calendar)

    # Plot the time series of significant wave height 
    fig, ax_hs = plt.subplots(figsize=(7.5,5))
    ax_hs.plot(awac4p5m_time, awac4p5m['waveHs'][:], color='r', label='4.5 m AWAC')
    ax_hs.plot(awac6m_time, awac6m['waveHs'][:], color='m', label='6 m AWAC')
    ax_hs.plot(array8m_time, array8m['waveHs'][:], color='y', label='8 m Array')
    
    # Set figure properties
    ax_hs.set_xlabel('Time [UTC]')
    ax_hs.set_ylabel('Significant Wave Height [m]')
    ax_hs.set_xlim(awac4p5m_time[0], awac4p5m_time[-1])
    ax_hs.set_ylim(0, 3.2)
    ax_hs.legend(loc=0)

    # Plot the time series of peak period
    
    # Add Mission Time block patches
    for ind in np.arange(1,len(start_times)):
        if ind == 6:
            # skip mission 6 which was not a real mission but a separate offload for the micros that were lost then recovered later - see notes spreadsheet
            continue
        ax_hs.add_patch(patches.Rectangle((start_times[ind], 0), end_times[ind]-start_times[ind], 3.5, linewidth=1, edgecolor='0.7', facecolor='0.7'))

    # Save the figure
    fig.savefig('../writing/Figures/fig02.png', dpi=300)

    return

if __name__=='__main__':
    main()