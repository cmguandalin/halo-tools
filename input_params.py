############################
# Parameters and constants #
############################

# Input cosmological parameters

h = 0.67556 
n_s = 0.9619
A_s = 2.215e-9
omega_cdm = 0.12038
omega_b = 0.022032
omega_k = 0.0
k_pivot = 0.05
N_ur = 3.046
N_ncdm = 0.0
T_cmb = 2.7255
w = -1.0

# Lazy parameters
pi = 3.14159265359
tpi = 2.*3.14159265359

# Derived parameters
Omega_m = (omega_b + omega_cdm)/(h**2.)
Omega_l = 1-Omega_m

# Speed of light
c = 299792.458 #km/s

# Gravitational constant in convenient units 
G = 4.302*pow(10,-9) #Mpc/Msun (km/s)^2

# Background critical density of the universe TODAY
rho_c = 3.0*pow(100,2.)/(8*pi*G) #h^2 Msun/Mpc^3 
# I have checked and it matches the PDG value 2.77e+11

# Background matter density TODAY
rho_m = Omega_m*rho_c 

# Still have to make this better (include on halo_model.py a function to compute M_star)
m_star = 1.5*1e12

# Threshold extrapolated in linear theory
delta_c = 1.686

# Local primordial non-Gaussianity amplitude @ CMB
fNL = 0.0
