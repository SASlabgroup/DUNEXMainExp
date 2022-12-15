"""
Tools to analyze individual waves from a mission
"""
import cftime
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from PyAstronomy import pyaC
from scipy import interpolate
from scipy import stats

def plot_mission_tracks(mission_dataset, bathy_file):
    """
    Plot all tracks from the mission over the bathymetry of the FRF

    Parameters
    ----------
    mission_dataset : nc.Dataset object
        The dataset object for the mission from a netCDF file
    bathy_file : str
        Path or url to the bathymetry file
    """
    fig, ax = plt.subplots(figsize=(8,8))

    # Create a bathymetry dataset from the bathymetry file and plot
    # a filled contour plot of the bathymetry
    bathy_dataset = nc.Dataset(bathy_file)
    xFRF_grid, yFRF_grid = np.meshgrid(bathy_dataset['xFRF'][:],
                                       bathy_dataset['yFRF'][:])
    bathy = bathy_dataset['elevation'][0,:,:]
    ax.contourf(xFRF_grid, yFRF_grid, bathy, cmap='gray')
    bathy_dataset.close()
    # cs = ax.contour(xFRF_grid, yFRF_grid, bathy, cmap='gray', linestyle=None)
    # ax.clabel(cs, levels=[-2,-4,-6], inline=True, fontsize=16, colors='k')

    # Plot the drift tracks of each microSWIFT on the mission
    x_locations = mission_dataset['xFRF']
    y_locations = mission_dataset['yFRF']
    time_vals = mission_dataset['time'][:]
    for n in range(mission_dataset['trajectory'].size):
        map = ax.scatter(x_locations[n,:], y_locations[n,:],
                         c=time_vals, cmap='plasma')

    # Setup the colorbar for the time axis
    map.set_clim([time_vals[0], time_vals[-1]])
    cbar = fig.colorbar(map, ax=ax, ticks=[time_vals[0], time_vals[-1]])
    time_labels = cftime.num2pydate([time_vals[0], time_vals[-1]],
                                    units=mission_dataset['time'].units,
                                    calendar=mission_dataset['time'].calendar)
    time_label_strs = [time_labels[0].strftime('%Y-%m-%d %H:%M'),
                       time_labels[1].strftime('%Y-%m-%d %H:%M')]
    cbar.ax.set_yticklabels(time_label_strs, rotation=0, va='center')
    cbar.ax.set_xlabel('Time [UTC]')

    # Figure properties
    ax.set_aspect('equal')
    ax.set_xlabel('Cross Shore Location [meters]')
    ax.set_ylabel('Along Shore Location [meters]')
    ax.plot([100,800],[505,505], linewidth=2, color='r')

def plot_mission_eta(mission_dataset):
    """
    Plot the time series of sea surface elevation for each microSWIFT
    on the mission.

    Parameters
    ----------
    mission_dataset : netCDF object
        netCDF4 object instance for the mission
    """
    fig, ax = plt.subplots()
    ax.set_xlabel('Time [UTC]')
    ax.set_ylabel('Sea Surface Elevation [m]')
    time = cftime.num2pydate(mission_dataset['time'],
                            units=mission_dataset['time'].units,
                            calendar=mission_dataset['time'].calendar)
    eta = mission_dataset['sea_surface_elevation']
    for n in range(mission_dataset['trajectory'].size):
        ax.plot(time, eta[n,:])

def compute_individual_waves(x_locations, y_locations, eta, time, bathy_file,
                             single_trajectory=False):
    """
    Compute a distribution of wave heights and their locations from
    Parameters
    ----------
    x_locations : np.ndarray
        Array of cross shore locations of microSWIFT
    y_locations : np.ndarray
        Array of along shore locations of microSWIFT
    eta : np.ndarray
        Array of sea surface elevation of microSWIFT
    time : np.ndaray
        Array of time values - not datetime values
    bathy_file : str
        path or url to the bathymetry file
    single_trajectory : boolean
        Whether or not we are plotting a single trajectory

    Returns
    -------
    wave_heights : list
        List of each wave height on a mission
    wave_x_locs : list
        List of each wave cross shore location
    wave_y_locs : list
        List of each wave along shore location
    """

    # Iniialize the arrays to store the wave height and locations
    wave_heights = []
    wave_x_locs = []
    wave_y_locs = []

    if single_trajectory is True:
        eta = eta.reshape(1,eta.size)
        x_locations = x_locations.reshape(1,x_locations.size)
        y_locations= y_locations.reshape(1,y_locations.size)
    else:
        pass

    for trajectory in range(eta.shape[0]):
        cross_time, cross_ind = pyaC.zerocross1d(time, eta[trajectory,:],
                                                getIndices=True)
        wave_inds = cross_ind[::2]

        # Compute Wave Height from each set of indices
        for n in np.arange(np.size(wave_inds)-1):
            # Get 45 elevation heights in between each zero crossing index
            eta_in_wave = eta[trajectory, wave_inds[n]:wave_inds[n+1]]
            x_in_wave = x_locations[trajectory, wave_inds[n]:wave_inds[n+1]]
            y_in_wave = y_locations[trajectory, wave_inds[n]:wave_inds[n+1]]

            # Set intital wave height to NaN, if there is a wave compute
            # the wave height and check if that overwrote the initial
            # NaN value, if it did then save the wave height and
            # locations
            wave_height = np.NaN
            if np.size(eta_in_wave) > 0:
                # Add the compute wave height
                wave_height = np.max(eta_in_wave) - np.min(eta_in_wave)
            else:
                pass

            if ~np.isnan(wave_height):
                wave_heights.append(wave_height)
                wave_x_locs.append(np.mean(x_in_wave))
                wave_y_locs.append(np.mean(y_in_wave))
            else:
                pass

    return wave_heights, wave_x_locs, wave_y_locs

def plot_wave_locations(wave_x_locs, wave_y_locs, bathy_file, color):
    """
    Plot the location of the individual waves from the mission

    Parameters
    ----------
    x_locations : list
        List of cross shore locations of waves
    y_locations : list
        List of along shore locations of waves
    bathy_file : list
        path or url to bathymetry file
    """
    fig, ax = plt.subplots(figsize=(8,8))

    # Create a bathymetry dataset from the bathymetry file and plot
    # a filled contour plot of the bathymetry
    bathy_dataset = nc.Dataset(bathy_file)
    xFRF_grid, yFRF_grid = np.meshgrid(bathy_dataset['xFRF'][:],
                                       bathy_dataset['yFRF'][:])
    bathy = bathy_dataset['elevation'][0,:,:]
    ax.contourf(xFRF_grid, yFRF_grid, bathy, cmap='gray')
    bathy_dataset.close()

    # Scatter the locations of the individual waves
    ax.scatter(wave_x_locs, wave_y_locs, color=color, s=1, marker='.')

    # Figure properties
    ax.set_aspect('equal')
    ax.set_xlabel('Cross Shore Location [meters]')
    ax.set_ylabel('Along Shore Location [meters]')
    ax.plot([100,800],[505,505], linewidth=2, color='r')

def plot_wave_height_dist(wave_heights, num_bins):
    """
    Plot the distibution of wave heights

    Parameters
    ----------
    wave_heights : list
        List of wave heights
    num_bins : int
        Number of bins for the histogram

    Returns
    -------
    ax : axis object
        axis object of the distribution
    """
    fig, ax = plt.subplots()

    # Fit a Raleigh distibution to the wave heights
    loc, scale = stats.rayleigh.fit(wave_heights)
    rayleigh_dist_x = np.linspace(min(wave_heights), max(wave_heights), 100)

    ax.hist(wave_heights, bins=num_bins, density=True,
            label=f'Number of Waves: {len(wave_heights)}')
    ax.plot(rayleigh_dist_x,
            stats.rayleigh(scale=scale, loc=loc).pdf(rayleigh_dist_x),
            label="Rayleigh Distribution Fit", color='r')
    ax.set_xlabel('Wave Height [m]')
    ax.set_ylabel('Probability Density [-]')
    ax.legend()

    return ax

def compute_sig_wave_height_var(wave_heights):
    """
    Compute Significant Wave Height, Hs from distribution of wave
    heights.

    Parameters
    ----------
    wave_heights : list
        List of wave heights

    Returns
    -------
    sig_wave_height : float
        significant wave height
    """
    sig_wave_height = 4 * np.sqrt(np.var(wave_heights))

    return sig_wave_height

def compute_sig_wave_height_top_third(wave_heights):
    """
    Compute Significant Wave Height, Hs from distribution of wave
    heights.

    Parameters
    ----------
    wave_heights : list
        List of wave heights

    Returns
    -------
    sig_wave_height : float
        significant wave height
    """
     # Sort wave heights
    wave_height_sort = np.sort(wave_heights)

    # Get Top 1/3rd of wave heights
    top_third_waves = wave_height_sort[-(len(wave_heights)//3):]

    # Average Top 1/3rd of wave heights
    sig_wave_height = np.mean(top_third_waves)

    return sig_wave_height

def compute_wave_bathy(x_locs, y_locs, bathy_file):
    """
    Compute the wave depth by interpolating the bathymetry data at
    the x and y locations of the individual waves.

    Parameters
    ----------
    x_locs : list
        List of cross shore locations of waves
    y_locs : list
        List of along shore locations of waves

    Returns
    -------
    wave_depths : list
        List of the depth of each wave
    """
    # Load in bathymetry Data and find location of each
    bathy = nc.Dataset(bathy_file)
    elevation = bathy['elevation'][0,:,:]
    bathy_xFRF = bathy['xFRF'][:]
    bathy_yFRF = bathy['yFRF'][:]
    bathy_f = interpolate.interp2d(bathy_xFRF, bathy_yFRF, elevation)

    # Depth of each wave
    wave_bathy = []
    for n in np.arange(len(x)):
        depth.append(np.squeeze(bathy_f(x[n], y[n])).item())
    return depth 