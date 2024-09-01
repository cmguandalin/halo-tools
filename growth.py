import numpy as np
from input_params import *

from scipy.interpolate import interp1d 
from scipy.integrate import odeint

# Define the Hubble parameter
# (NEGLECTING radiation and assuming that dark energy is a cosmological constant (w=-1)!)
def Ha(a,h,om,ol):
    return 100*h*np.sqrt(om*a**(-3.) + ol + (1-ol-om)*a**(-2.))

# Define the derivative of the Hubble parameter wrt a:
def derH(a,h,om,ol):
    return (50*h*(-((2*(1-ol-om))/a**3.)-(3*om)/a**4.))/np.sqrt(ol+(1-ol-om)/a**2.+om/a**3.)

# Define the function that yields the system of two frist order equations we would like to solve
def g(v,a,h,om,ol):
    delta, y = v
    dvdt = [y, -(derH(a,h,om,ol)/Ha(a,h,om,ol) + 3/a)*y + (1.5*om/(a**5.))*((Ha(a,h,om,ol)/(100*h))**(-2.))*delta]
    return dvdt

# Initial conditions: v0 = [delta(a_in), y=delta'(a_in)]
# We set delta = a at early (matter dominated) times, so delta' = 1.
v0 = [0.00000001,1.]

# Time steps: an array of values starting from 0.01 going up to 1., in steps of 0.01.
# With this, we will have at hands the solution with the scale factor from 0.01 to 1.
a = np.arange(0.00000001, 1.1, 0.00001)

def D(z):
    # Solve the equation with the help of odeint solver
    delta_lcdm = odeint(g, v0, a, (h, Omega_m, Omega_l))

    z_plot = 1/a - 1.
    d_lcdm = interp1d(z_plot,delta_lcdm[:,0])

    return d_lcdm(z)/d_lcdm(0.0)
