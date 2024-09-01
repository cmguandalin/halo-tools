import numpy as np
from input_params import *
import distances as dist

####################
# Useful functions #
####################

# Obs: rho_m is the background matter density today (defined at input_params).
#      mass_to_radius thus returns a comoving scale.
#      with these definitions, the dn/dlnM will be given in terms of comoving volume number counts.

def E_sq(z,om,w):
    return pow(dist.H(z,om,w)/100,2.0)

# Background matter density at redshift z
def rho_mz(z):
    return rho_m*(1+z)**3.

def rho_cz(z,om,w):
    return rho_c*E_sq(z,om,w)

def omega_m(z,om,w):
    return om*pow(1+z,3.0)/E_sq(z,om,w)

def Delta_c(z,om,w):
    x = omega_m(z,om,w) - 1.0
    return 18*(pi**2.0) + 82*x - 39*(x**2.0)

def Delta(z,om,w):
    return Delta_c(z,om,w)/om

def r_vir(z,om,w,M):
    return pow(((3./4.)*(M/pi))/(rho_cz(z,om,w)*Delta_c(z,om,w)),1./3.)

def c_fun(z,M,Mstar):
    return (9.0/(1+z))*pow(M/Mstar,-0.13)

def rho_s(z,om,w,M,Mstar):
    rvir = r_vir(z,om,w,M)
    conc = c_fun(z,M,Mstar)
    factor = 4*pi*pow(rvir/conc,3.)
    term = np.log(1+conc)-(conc/(1+conc))
    return (M/factor)/term

def mass_to_radius(M):
    return pow(3*M/(4*pi*rho_m),1./3.) #Mpc/h

def radius_to_mass(R):
    return 4.0*pi*rho_m*R**3/3.0 

def win_fun(k,R):
    return 3.0*(np.sinc(k*R/pi)-np.cos(k*R))/pow(k*R,2.0)

def integrand_win_fun(k,R):
    return (pow(k*win_fun(k,R),2.))/(2*pow(pi,2.))
