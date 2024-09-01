import numpy as np

############################
# Galaxy Cluster Abundance #
############################


'''
    Multiplicity functions only

'''

# Tinker fits at ANY z:
def alpha_fun(D):
    log_alpha = -pow(0.75/(np.log10(D/75)),1.2)
    return pow(10.0,log_alpha)

def tinker2008_f(z,D):
    # Redshift dependence
    alpha = alpha_fun(D)
    z_A = (1+z)**(-0.14)
    z_a = (1+z)**(-0.06)
    z_b = (1+z)**(-alpha)

    # Parameters A, a, b, c

    if D == 200:
        params = [0.186*z_A, 1.47*z_a, 2.57*z_b, 1.19]
    elif D == 300:
        params = [0.200*z_A, 1.52*z_a, 2.25*z_b, 1.27]
    elif D == 400:
        params = [0.212*z_A, 1.56*z_a, 2.05*z_b, 1.34]
    elif D == 600:
        params = [0.218*z_A, 1.61*z_a, 1.87*z_b, 1.45]
    elif D == 800:
        params = [0.248*z_A, 1.87*z_a, 1.59*z_b, 1.58]
    elif D == 1200:
        params = [0.255*z_A, 2.13*z_a, 1.51*z_b, 1.80]
    elif D == 1600:
        params = [0.260*z_A, 2.30*z_a, 1.46*z_b, 1.97]
    elif D == 2400:
        params = [0.260*z_A, 2.53*z_a, 1.44*z_b, 2.24]
    elif D == 3200:
        params = [0.260*z_A, 2.66*z_a, 1.41*z_b, 2.44]
    else:
        params = [0, 0, 0, 0, 0]
        print('Enter a valid Delta. Your value:', D)
    
    return np.array(params)

def tinker2008_g(D):
    # Parameters B, d, e, f, g
    if D == 200:
        params = [0.482, 1.97, 1.00, 0.51, 1.228]
    elif D == 300:
        params = [0.466, 2.06, 0.99, 0.48, 1.310]
    elif D == 400:
        params = [0.494, 2.30, 0.93, 0.48, 1.403]
    elif D == 600:
        params = [0.494, 2.56, 0.93, 0.45, 1.553]
    elif D == 800:
        params = [0.496, 2.83, 0.96, 0.44, 1.702]
    elif D == 1200:
        params = [0.450, 2.92, 1.04, 0.40, 1.907]
    elif D == 1600:
        params = [0.466, 3.29, 1.07, 0.40, 2.138]
    elif D == 2400:
        params = [0.429, 3.37, 1.12, 0.36, 2.394]
    elif D == 3200:
        params = [0.388, 3.30, 1.16, 0.33, 2.572]
    else:
        params = [0, 0, 0, 0, 0]
        print('Enter a valid Delta. Your value:', D)
    
    return np.array(params)

# Multiplicity functions

def tinker_f(z,D,s):
    A=tinker2008_f(z,D)[0]
    a=tinker2008_f(z,D)[1]
    b=tinker2008_f(z,D)[2]
    c=tinker2008_f(z,D)[3]
    
    return A*(pow(s/b,-a)+1)*np.exp(-c/pow(s,2.0))

def tinker_g(D,s):
    B=tinker2008_g(D)[0]
    d=tinker2008_g(D)[1]
    e=tinker2008_g(D)[2]
    f=tinker2008_g(D)[3]
    g=tinker2008_g(D)[4]
    
    return B*(pow(s/e,-d)+pow(s,-f))*np.exp(-g/pow(s,2.0))

def PS(s):
    nu = 1.686/s
    
    return np.sqrt(2.0/np.pi)*nu*np.exp(-pow(nu,2.0)/2.0)

def ST(s):
    A = 0.322
    q = 0.707
    p = 0.300
    
    nu = 1.686/s
    
    return A*np.sqrt(2.0*q/np.pi)*(1+pow(q*pow(nu,2.0)/2.0,-p))*nu*np.exp(-q*pow(nu,2.0)/2.0)
    
