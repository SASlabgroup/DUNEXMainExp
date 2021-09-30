# Coastal Engineering Python Module File 
# Written By: EJ Rainville, Spring 2021
# Description: This script contains functions and defintions from the course CEWA573. 

# Read in Functions 
# Define read CDIP function
def readCDIP(filename):
    file = open(filename, 'r')
    lines = file.readlines()

    # Loop through file and split each row 
    lines_split = []
    for n in np.arange(len(lines)):
        lines_split.append(lines[n].split())

    # Organize the data into vectors for each variable
    # Dates for data
    dates = pd.to_datetime([''.join(item[:5]) for item in lines_split[2:]])

    # Significant Wave Height, units are meters
    Hs = [float(item[5]) for item in lines_split[2:]] 

    # Peak Period, units are seconds
    Tp = [float(item[6]) for item in lines_split[2:]]

    # Peak Direction, units are degrees
    Dp = [float(item[7]) for item in lines_split[2:]]
    
    # Return Variables
    return dates, Hs, Tp, Dp
    

# Compute Wavenumber from Dispersion Relation 
def wavenumber(f, depth, g=9.81):
    import numpy as np
    """ input frequency (Hz) and depth (m)
        converted from Dr. Jim Thomsons wavelength.m
        adapted with another student as well
    """
    omega = 2*np.pi*f 
    depthR = np.round(depth)
    if depth<20: #shallow
        guess_k = np.sqrt(omega/(g*depthR))
        eps = 0.01*guess_k 
        err = np.abs(omega**2 - g*guess_k*np.tanh(guess_k*depthR))
    else:
        guess_k = omega**2/g
        eps = 0.01*guess_k
        err = np.abs(omega**2 - g*guess_k*np.tanh(guess_k*depthR))
    k = guess_k
    while err>eps:
        k = guess_k - (omega**2 - g*guess_k*np.tanh(guess_k*depthR))/(-g*np.tanh(guess_k*depthR) - g*guess_k*depthR*np.cosh(guess_k)**2)
        err = np.abs(omega**2 - g*k*np.tanh(k*depthR))
        guess_k = k
    return k

# Ray Tracing Functions
# Written By: Emma Nuss, Adapted from Jim Thomson
def ray_tracing(f, alpha0, depths, xi, x, y):
    """ Ray tracing function for alongshore uniform bathymetry
        Takes inputs of frequency (f) in Hz, inital angle (alpha0) in radians,
        cross-shore profile of depths (depths), x points (xi) where np.max(xi) 
        gives shoreline x position, and initial (x,y) coordinates
        Output:
        x and y points of ray and alpha at each time-step
    """
    omega = 2*np.pi*f 
    k = np.asarray([wf.wavenumber(f, d) for d in depths])
    c = np.asarray([wf.phase_speed(omega, k[i], depths[i]) for i in range(len(depths))])
    dx = 1
    alpha = alpha0
    alpha_all = []
    alpha_all.append(alpha)
    x_all = []
    x_all.append(x)
    y_all = []
    y_all.append(y)

    dt = 0.1
    while x < np.max(xi):
        ct = np.interp(x,yi,c)
        alpha = np.arcsin(ct/c[0]*np.sin(alpha0))
        x = x + ct*dt*np.cos(alpha)
        y = y + ct*dt*np.sin(alpha)
        alpha = np.arcsin(ct/c[0]*np.sin(alpha0))
        alpha_all.append(alpha) 
        x_all.append(x)
        y_all.append(y)
    return np.asarray(x_all), np.asarray(y_all), np.asarray(alpha_all)*180/np.pi
