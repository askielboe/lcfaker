# -*- coding: utf-8 -*-
"""
File:    synthesizers.py

Author:  Andreas Skielboe (skielboe@dark-cosmology.dk)
Date:    March 2013

Summary:

    Functions for generating lightcurves, using data or synthesizers.

"""

import numpy as np
from LightCurve import LightCurve

def getLightCurveFromSpectra(spectra, minimum, maximum):
    """
    Generates the lightcurve of integrated flux within limits defined by minimum and maximum
    For continuum, could be (from previous getLightCurveCont() function):
    minimum = 5183
    maximum = 5193
    """
    time = np.zeros(len(spectra))
    flux = np.zeros(len(spectra))
    ferr = np.zeros(len(spectra))

    for i,spectrum in enumerate(spectra):
        time[i] = spectrum.date
        flux[i], ferr[i] = spectrum.integrate(minimum, maximum)

    return LightCurve(time, flux, ferr)

def getLightCurveFromSpectraHBeta(spectra):
    """
    Generates a lightcurve from integrating the HBeta line in a series of spectra.
    Input: List of spectra instances.
    Output: Lightcurve instance.
    """
    time = np.zeros(len(spectra))
    flux = np.zeros(len(spectra))
    ferr = np.zeros(len(spectra))

    for i, spectrum in enumerate(spectra):
        time[i] = spectrum.date
        flux[i], ferr[i] = spectrum.integrateHBeta()

    return LightCurve(time, flux, ferr)

def syntheticLightCurveWiener(nBins):
    """
    Function to generate random walk light curve using pyprocesses.py
    """
    from lcfaker.vendor import pyprocess
    from numpy import asarray

    wienerParams = {"mu":0, "sigma":0.1}
    wienerInitial = {"startTime":0, "startPosition":0, "endTime":nBins, "endPosition": 0 }
    Wiener = pyprocess.Wiener_process(wienerParams, wienerInitial)

    lc = asarray(Wiener.generate_sample_path(range(nBins)))

    return LightCurve(lc.T[0], lc.T[1])

def syntheticLightCurveOU(nBins):
    """
    Function to generate random walk light curve using pyprocesses.py
    """
    from lcfaker.vendor import pyprocess
    from numpy import asarray

    OUParams = {"theta":1., "mu":5., "sigma":1.25}
    OUInitial = {"startTime":0, "startPosition":0, "endTime":nBins, "endPosition": 0 }
    OU = pyprocess.OU_process(OUParams, OUInitial)

    lc = asarray(OU.generate_sample_path(range(nBins)))

    return LightCurve(lc.T[0], lc.T[1])

def syntheticLightCurveMacLeod(nDays = 200, maxDays = 200.0, \
                                magApparant = 19, magAbsolute = -24.0, \
                                mass = 1.0e9, \
                                lambdarf = 5100.0, \
                                z = 0.5, dist = 1.0e6):
    """
    Function to generate lightcurve based on damped random walk
    using equations from MacLeod et al.
    """
    import numpy as np
    from lib.ouprocess import OUprocess
    from lib.units import magi_to_fluxJy

    print "Generating synthetic light curve using method described in MacLeod et al."
    print "Parameters:"
    print "nDays = ", nDays
    print "maxDays = ", maxDays
    print "magApparant = ", magApparant
    print "magAbsolute = ", magAbsolute
    print "mass = ", mass
    print "lamdarf = ", lambdarf
    print "z = ", z
    print "dist = ", dist

    magApparant = 5.0 * np.log10(dist) - 5.0 + magAbsolute

    ### Calculate random walk parameters from physical parameters
    # Relation from MacLeod et al. 2010 eq. 7
    # Parameters from MacLeod et al. 2010 Table 1
    A = -0.51
    B = -0.479
    C = 0.131
    D = 0.18
    E = 0.0

    logsf = A \
        + B * np.log(lambdarf/4000.) \
        + C * (magAbsolute + 23.) \
        + D * np.log(mass/1e9) \
        + E * np.log(1 + z)

    sf = np.exp(logsf)

    A = 2.4
    B = 0.17
    C = 0.03
    D = 0.21
    E = 0.0

    logtau = A \
        + B * np.log(lambdarf/4000.) \
        + C * (magAbsolute + 23.) \
        + D * np.log(mass/1e9) \
        + E * np.log(1 + z)

    tau = np.exp(logtau)
    ###

    time, mag = OUprocess(nDays, maxDays, magApparant, sf, tau)

    # Convert from magnitude to flux
    flux = magi_to_fluxJy(mag)

    return LightCurve(time, flux)
