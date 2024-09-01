import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import simps

from input_params import *
import utilities as utils
import growth as growth
import multiplicity as mult



class halo_tools(object):
    """

    	This class computes the following halo properties:
		- halo mass function
		- halo bias
	    
        k: vector of wavenumbers corresponding the P(k) (units of h/Mpc)
        P: linear matter power-spectrum at z=0 (units of Mpc^3/h^3)

	Obs: if you want to change Delta = 200, change directly in lines
	     Ln 59 and Ln 62, and for the bias, Ln 82
        
    
    """
    def __init__(self, k, P):

        self.k = k
        self.P = P
        self.Plin = interp1d(self.k,self.P)

        self.r_min = 2.*np.pi/max(k)
        self.r_max = 2.*np.pi/min(k)
        # print 'The minimum scale (corresponding to k_max) we can go is: R_min =', r_min
        self.r_int = np.logspace(np.log10(self.r_min),np.log10(self.r_max),10000)
        self.mass_int = utils.radius_to_mass(self.r_int)
        # print 'The minimum mass (corresponding to r_min) we have is: M_min =', min(mass_int)
    
        self.var = []
        for i in self.r_int:
            self.var.append(simps(utils.integrand_win_fun(self.k,i)*self.Plin(self.k),self.k))
        
        self.var_tmp = np.array(self.var)

        self.variance = interp1d(self.r_int,np.sqrt(self.var_tmp),bounds_error=False,fill_value=0.0)
        self.dlns_dlnm = interp1d(self.mass_int,-(self.mass_int/self.variance(self.r_int))*np.gradient(self.variance(self.r_int),self.mass_int))



    def mass_function(self, M, z, f):
        self.M = M
        self.z = z
        self.R = utils.mass_to_radius(self.M)
        self.sigma = growth.D(self.z)*self.variance(self.R)
        self.ln_sigma = np.log(1.0/self.sigma)

        if f == 0:
            self.f_sigma = mult.tinker_f(self.z,200,self.sigma)
    
        elif f == 1:
            self.f_sigma = mult.tinker_g(200,self.sigma)
        
        elif f == 2: 
            self.f_sigma = mult.PS(self.sigma)
        
        else:     
            self.f_sigma = mult.ST(self.sigma)
        
        return self.f_sigma*(rho_m/self.M)*self.dlns_dlnm(self.M)


    def bias(self, M, z, f):
        self.M = M
        self.z = z
        self.R = utils.mass_to_radius(self.M) 
        self.NU = 1.686/(growth.D(self.z)*self.variance(self.R))

        if f == 0 or f == 1:
            # Parameters
            self.y = np.log10(200)
            self.A = 1.0 + 0.24*self.y*np.exp(-pow(4/self.y,4.0))
            self.a = 0.44*self.y - 0.88
            self.B = 0.183
            self.b = 1.500
            self.C = 0.019 + 0.107*self.y + 0.19*np.exp(-pow(4/self.y,4.0))
            self.c = 2.4
            self.resp = 1.0 - self.A*pow(self.NU,self.a)/(pow(self.NU,self.a)+(1.686**self.a)) + self.B*pow(self.NU,self.b) + self.C*pow(self.NU,self.c)
    
        elif f == 2:
            self.resp = 1.0 + (pow(self.NU,2.0)-1.0)/1.686
        
        elif f == 3:
            self.a = 0.707
            self.p = 0.3
            self.resp = 1.0 + (self.a*pow(self.NU,2.0) - 1.0)/1.686 + 2.0*self.p/(1.686*(1+pow(self.a*(self.NU**2.0),self.p)))
    
        return self.resp