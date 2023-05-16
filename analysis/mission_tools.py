"""
Tools to analyze individual waves from a mission
"""
import cftime
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from PyAstronomy import pyaC
from scipy import interpolate
from scipy import signal
from scipy import stats

def plot_mission_tracks(mission_dataset, bathy_file, trajectory_subset=None,
                        time_ind_subset=None):
    """
    Plot all tracks from the mission over the bathymetry of the FRF

    Parameters
    ----------
    mission_dataset : nc.Dataset object
        The dataset object for the mission from a netCDF file
    bathy_file : str
        Path or url to the bathymetry file
    trajectory_subset : list
        A list of trajectory indices that you want to plot.
    trajectory_subset : np.ndarray
        An array of time indices that you want to plot.
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
    if trajectory_subset is None and time_ind_subset is None:
        x_locations = mission_dataset['xFRF']
        y_locations = mission_dataset['yFRF']
        time_vals = mission_dataset['time'][:]
        for n in range(mission_dataset['trajectory'].size):
            map = ax.scatter(x_locations[n,:], y_locations[n,:],
                            c=time_vals, cmap='plasma')

    elif trajectory_subset is not None and time_ind_subset is None:
        x_locations = mission_dataset['xFRF']
        y_locations = mission_dataset['yFRF']
        time_vals = mission_dataset['time'][:]
        for n in range(mission_dataset['trajectory'].size):
            if n in trajectory_subset:
                map = ax.scatter(x_locations[n,:], y_locations[n,:],
                            c=time_vals, cmap='plasma')
            else:
                pass

    if trajectory_subset is None and time_ind_subset is not None:
        x_locations = mission_dataset['xFRF']
        y_locations = mission_dataset['yFRF']
        time_vals = mission_dataset['time'][:]
        time_ind_index = 0
        for n in range(mission_dataset['trajectory'].size):
            map = ax.scatter(x_locations[n,time_ind_subset[time_ind_index]],
                                y_locations[n,time_ind_subset[time_ind_index]],
                                c=time_vals[time_ind_index], cmap='plasma')
            time_ind_index += 1

    elif trajectory_subset is not None and time_ind_subset is not None:
        x_locations = mission_dataset['xFRF']
        y_locations = mission_dataset['yFRF']
        time_vals = mission_dataset['time'][:]
        time_ind_index = 0
        for n in range(mission_dataset['trajectory'].size):
            if n in trajectory_subset:
                map = ax.scatter(x_locations[n,time_ind_subset[time_ind_index]],
                                 y_locations[n,time_ind_subset[time_ind_index]],
                                 c=time_vals[time_ind_subset[time_ind_index]],
                                 cmap='plasma')
                time_ind_index += 1
            else:
                pass

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
    awac4p5m_location = [397.35, 890.98] # Converted from lat lon locations
                                         # published on FRF data portal
    ax.scatter(awac4p5m_location[0],awac4p5m_location[1],
               color='r', label='4.5 m AWAC')
    ax.set_aspect('equal')
    ax.set_xlabel('Cross Shore Location [meters]')
    ax.set_ylabel('Along Shore Location [meters]')
    ax.plot([50,591],[510,510], linewidth=2, color='r', label='Pier')
    ax.legend()
    
    return ax

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

    return ax

def compute_individual_waves(x_locations, y_locations, eta, time, bathy_file,
                             single_trajectory=False, time_step=1/12):
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
    time_step : float
        Value for the difference in time between two points

    Returns
    -------
    wave_heights : list
        List of each wave height on a mission
    wave_periods : list
        List of each wave period on a mission
    wave_x_locs : list
        List of each wave cross shore location
    wave_y_locs : list
        List of each wave along shore location
    """

    # Initialize the arrays to store the wave height and locations
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
            # Get elevation heights in between each zero crossing index
            eta_in_wave = eta[trajectory, wave_inds[n]:wave_inds[n+1]]
            x_in_wave = x_locations[trajectory, wave_inds[n]:wave_inds[n+1]]
            y_in_wave = y_locations[trajectory, wave_inds[n]:wave_inds[n+1]]

            # Set initial wave height to NaN, if there is a wave compute
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
    ax.plot([50,591],[510,510], linewidth=2, color='r')

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

def compute_sig_wave_height_var(eta, single_trajectory=False):
    """
    Compute Significant Wave Height, Hs from time series of sea surface
    elevation.

    Parameters
    ----------
    eta : np.ndarray
        2D array of sea surface elevation
    single_trajectory : boolean
        True or False if this is a single trajectory

    Returns
    -------
    sig_wave_height : float
        significant wave height
    """
    if single_trajectory is True:
        eta = eta.reshape(1,eta.size)
    else:
        pass

    # compute varaince for each trajectory and average across each traj
    variance_on_each_trajectory = np.nanvar(eta, axis=1)
    average_variance = np.mean(variance_on_each_trajectory)
    sig_wave_height = 4 * np.sqrt(average_variance)

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
    standard_dev : float
        standard deviation of wave heights
    """
     # Sort wave heights
    wave_height_sort = np.sort(wave_heights)

    # Get Top 1/3rd of wave heights
    top_third_waves = wave_height_sort[-(len(wave_heights)//3):]

    # Average Top 1/3rd of wave heights
    sig_wave_height = np.mean(top_third_waves)

    standard_dev = np.std(top_third_waves)

    return sig_wave_height, standard_dev

def compute_sig_wave_height_rms(wave_heights):
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
    h_rms = np.sqrt(np.mean(np.square(wave_heights)))
    sig_wave_height = 1.416 * h_rms

    return sig_wave_height

def closest_awac_sig_wave_height(mission_time, awac_file):
    """
    Find the closest AWAC significant wave height

    Parameters
    ----------
    mission_time : float
        time of the mission you are comparing to the awac
    awac_file : str
        path or url to the awac file

    Returns
    -------
    awac_sig_wave_height : float
        significant wave height
    """
    awac_data = nc.Dataset(awac_file)

    awac_sig_wave_height = np.interp(mission_time,
                        awac_data['time'][:],
                        awac_data['waveHs'][:])

    awac_data.close()
    return awac_sig_wave_height

def closest_awac_spectra(mission_time, awac_file):
    """
    Find the closest AWAC spectrum

    Parameters
    ----------
    mission_time : float
        time of the mission you are comparing to the awac
    awac_file : str
        path or url to the awac file

    Returns
    -------
    awac_spectra : np.ndarray
        energy density spectrum
    awac_freq : np.ndarray
        frequencies associated with the spectrum
    """
    awac_data = nc.Dataset(awac_file)

    awac_freq = awac_data['waveFrequency'][:]

    awac_spectra = np.empty(awac_data['waveEnergyDensity'][:].shape[1])
    
    for n in range(awac_spectra.size):
        awac_spectra[n] = np.interp(mission_time,
                                    awac_data['time'][:],
                                    awac_data['waveEnergyDensity'][:,n])

    awac_data.close()
    return awac_spectra, awac_freq

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
    for n in np.arange(len(x_locs)):
        wave_bathy.append(np.squeeze(bathy_f(x_locs[n], y_locs[n])).item())
    return wave_bathy

def bathy_along_track(bathy_file:str, xFRF:np.ndarray, yFRF:np.ndarray,
                      single_trajectory=False):
    """
    Linearly interpolates the bathymetry along the track of
    the microSWIFT.

    Parameters
    ----------
    bathy_file : str
        url or path to bathy bathymetry file
    xFRF : np.ndarray
        1D or 2D array of microSWIFT xFRF locations
    yFRF : np.ndarray
        1D or 2D array of microSWIFT xFRF locations
    single_trajectory : boolean
        True or False if plotting a single trajectory

    Returns
    -------
    bathy_along_track : np.ndarray
        1D or 2D array of bottom elevation at each location along the track

    """
    if single_trajectory is True:
        xFRF = xFRF.reshape(1,xFRF.size)
        yFRF = yFRF.reshape(1,yFRF.size)
    else:
        pass
    
    # Create bathymetry interpolating function from 2D grid
    bathy_dataset = nc.Dataset(bathy_file)
    bathy_xFRF = bathy_dataset['xFRF'][:]
    bathy_yFRF = bathy_dataset['yFRF'][:]
    bathy = bathy_dataset['elevation'][0,:,:]
    bathy_f = interpolate.interp2d(bathy_xFRF, bathy_yFRF, bathy)

    bathy_along_track = np.empty(xFRF.shape)
    for trajectory in range(xFRF.shape[0]):
        for n in np.arange(xFRF.shape[1]):
            bathy_along_track[trajectory, n] = np.squeeze(
                                                bathy_f(xFRF[trajectory, n],
                                                yFRF[trajectory, n]).item())

    return np.array(bathy_along_track)

def ind_in_depth(bathy_along_tracks, depth_min, depth_max,
                 single_trajectory=False):
    """
    Return the indices that the depth along the track is between the
    minimum and maximum values given.

    Parameters
    ----------
    bathy_along_track : np.ndarray
        Bathymetry along the tracks
    depth_min : float
        minimum depth you are interested in, minimum is the most
        negative value
    depth_max : float
        maximum depth you are interested in, maximum is the most
        positive value
    single_trajectory : boolean
        True or False if plotting a single trajectory

    Returns
    -------
    in_depth_indices : list
        list of largest arrays of indices in the depth range
    """
    if single_trajectory is True:
        bathy_along_tracks = bathy_along_tracks.reshape(1,
                                                       bathy_along_tracks.size)
    else:
        pass

    in_depth_indices = []

    for n in np.arange(bathy_along_tracks.shape[0]):
        bathy = bathy_along_tracks[n,:]
        inds_in_set = np.where((bathy > depth_min)
                               & (bathy < depth_max))[0]

        arrays = np.split(inds_in_set, np.where(np.diff(inds_in_set) != 1)[0]+1)
        size = 0
        ind_of_max = 0
        for n in range(len(arrays)):
            if len(arrays[n]) > size:
                size = len(arrays[n])
                ind_of_max = n

        in_depth_indices.append(arrays[ind_of_max])

    return in_depth_indices

def compute_spectra(z:np.ndarray, fs:float):
    """
    Compute energy density spectrum from sea surface elevation
    time series.

    Parameters
    ----------
    z : np.ndarray
        Array that contains time series you want to compute the energy
        spectrum for.
    fs : float
        sampling frequency for the time series
    """
    nperseg = 3600 # Window size in points
    overlap = 0.50
    f_raw, E_raw = signal.welch(z, fs=fs, window='hann', nperseg=nperseg,
                                noverlap=np.floor(nperseg*overlap))

    # Band Average the Spectra
    points_to_average = 5
    num_sections = E_raw.size // points_to_average
    f_chunks = np.array_split(f_raw, num_sections)
    E_chunks = np.array_split(E_raw, num_sections)
    f_total = np.array([np.mean(chunk) for chunk in f_chunks ])
    E_total = np.array([np.mean(chunk) for chunk in E_chunks ])
    dof = np.floor(points_to_average * (8/3) * (z.shape[0]/ (nperseg//2)))

    inds_to_return = np.where((f_total > 0.02) & (f_total < 0.5))
    f = f_total[inds_to_return]
    E = E_total[inds_to_return]

    return f, E, dof
