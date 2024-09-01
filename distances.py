import numpy as np
from scipy.interpolate import interp1d 
from scipy import integrate
from scipy import interpolate
from scipy import special

from input_params import *

##################################
# Functions to compute distances #
##################################

# Hubble parameter as a function of redshift, in units of (h/Mpc) (Km/s)
def H(z,om,w):
    return 100*np.sqrt(om*(1+z)**3. + (1-om)*(1+z)**(3.*(1+w))) 

# Comoving distance, in units of Mpc/h
def comov(z,om,w):
    
    if isinstance(z, (list, tuple, np.ndarray)) == True:
        comov = np.zeros(len(z))
        integrand = lambda z: 1./H(z,om,w)
        for i in range(len(z)):
            integral = integrate.quad(integrand,0.0,z[i])
            comov[i] = c*integral[0]
    else:
        integrand = lambda z: 1./H(z,om,w)
        integral = integrate.quad(integrand,0.0,z)
        comov = c*integral[0]
        
    return comov

# Luminosity distance, in units of Mpc/h
def lum(z,om,w):
    
    if isinstance(z, (list, tuple, np.ndarray)) == True:
        lum_dist = np.zeros(len(z))
        integrand = lambda z: 1./H(z,om,w)
        for i in range(len(z)):
            integral = integrate.quad(integrand,0.0,z[i])
            lum_dist[i] = c*integral[0]*(1+z[i])
    else:
        integrand = lambda z: 1./H(z,om,w)
        integral = integrate.quad(integrand,0.0,z)
        lum_dist = c*integral[0]*(1+z)
        
    return lum_dist

# Angular distance (equals to the comoving distance), in units of Mpc/h
# We keep this function since we define some things depending on it. Must clean this up later...
def ang(z,om,w):
    ang_dist = np.copy(comov(z,om,w))
    return ang_dist

# Physical angular distance, in units of Mpc/h
def phys_ang(z,om,w):
    
    if isinstance(z, (list, tuple, np.ndarray)) == True:    
        phys_ang_dist = np.zeros(len(z))
        r_obs = comov(z,om,w)
        for i in range(len(z)):
            phys_ang_dist[i] = r_obs[i]/(1+z[i])
            
    else: 
        phys_ang_dist = comov(z,om,w)/(1+z)
        
    return phys_ang_dist
