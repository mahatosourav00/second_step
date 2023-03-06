# do only hapke fit and save derived parameter

import numpy as np
from spectral import *
from scipy.optimize import curve_fit
import os



def load_flux(fluxname):

    with open(fluxname, 'r') as file:
        data = file.readlines()
        size = len(data)
        flux = np.zeros(size)
        wavelengths = np.zeros(size)
        for i in range(size):
            flux[i]=(float(data[i].strip().split(" ,")[1]))
            wavelengths[i]=(float(data[i].strip().split(",")[0]))
    return flux,wavelengths

def single_term_HG(alpha, g):
    return ((1-g**2)/(1+(2*g*np.cos(np.radians(alpha)))+g**2)**(3/2))

def double_term_HG(alpha, b, c):
    return ((((1+c)/2)*single_term_HG(alpha, -b)) + (((1-c)/2)*single_term_HG(alpha, b)))

#shadow hiding opposition effect function----
def shoe(alpha, Bs0, hs):
    return (Bs0/(1+((1/hs)*np.tan(np.radians(alpha/2)))))

#isotropic multiple scattering
def chandrasekhar_H(w,x):
    gamma = np.sqrt(1-w)
    r0 = (1-gamma)/(1+gamma)
    sec_part = w*x*(r0 + ((1-2*r0*x)/2) * np.log((1+x)/x))
    return ((1 - sec_part)**(-1))   

#hapke modelling------
#with 2-HG phase funtion-----
def hapke_single_scattering(X, w, b, c, Bs0, hs): # returns reflected brightness
    #define geometry
    i, e, alpha = X
    #define phase function
    p = double_term_HG(alpha, b, c)
    #define shadow hiding effect
    Bsh = shoe(alpha, Bs0, hs)
    #p = cornette_shanks(alpha, g)
    mu0 = np.cos(np.radians(i))
    mu = np.cos(np.radians(e))

    H0 = chandrasekhar_H(w,mu0)  
    H = chandrasekhar_H(w,mu)

    return ((w/(4*np.pi))*(mu0/(mu0+mu))*(p*(1 + Bsh)) + ((H0 * H) - 1)) #It is only for 40th band

def Hapke_fit(i, e, alpha, obs_ref):
    popt, pcov = curve_fit(hapke_single_scattering,(i, e, alpha), obs_ref, bounds = [(-np.inf, 0, -2, 0, 0),(1, 1, 2, np.inf, np.inf)])
    return popt, pcov


def main(binned, band):

    flux_dir = 'm3_solarflux.txt'
    flux, wavelengths = load_flux(flux_dir)
    distance = 0.983750743393
    flux = flux/(np.pi*(distance)**2)
    wavelengths = np.array(wavelengths[:85])/1e3

    alpha = []
    e = []
    i = []
    obs_ref = []
    for ph in range(len(binned)):
        for em in range(len(binned[0])):
            for inc in range(len(binned[0,0])):
                if binned[ph, em, inc][1] == 0:
                    pass
                else:
                    alpha.append(ph)
                    e.append(em)
                    i.append(inc)
                    obs_ref.append(binned[ph,em,inc][0]/binned[ph,em,inc][1])
    alpha = np.array(alpha)
    e = np.array(e)
    i = np.array(i)
    obs_ref = np.array(obs_ref)/flux[band]
    popt, pcov = Hapke_fit(i, e, alpha, obs_ref)
    result = (popt, pcov)
    if not os.path.exists('Hapke_parameters'):
        os.makedirs('Hapke_parameters')
    np.save(os.path.join('Hapke_parameters', 'params_%d.npy'%band), result)

