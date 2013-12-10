# -*- coding: utf-8 -*-
"""
File:    dataLoaders.py

Author:  Andreas Skielboe (skielboe@dark-cosmology.dk)
Date:    May 2013

Summary:

    Data loaders.

"""
import numpy as np
from Spectrum import Spectrum

def loadSpectrumWandersPeterson(filename):
    name = filename.split('/')[-1]
    date = name[2:6]
    data = np.loadtxt(filename).T
    spectrum = Spectrum(wavelength = data[0], flux = data[1], ferr = data[2], date = date)
    spectrum.lineWindow = (4870, 5010)
    spectrum.contWindow = [(4600,4740),(4790,4840),(5130,5300)]
    return spectrum

def loadSpectrumSDSS(filename):
	wavelength, flux, bestfit, skyflux = np.loadtxt(filename, delimiter=',', skiprows=1).T
	spectrum = Spectrum(wavelength = wavelength, flux = flux)
	return spectrum
