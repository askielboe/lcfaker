# -*- coding: utf-8 -*-
"""
File:    reverberation.py

Author:  Andreas Skielboe (skielboe@dark-cosmology.dk)
Date:    March 2013

Summary:

    Functions for doing reverberation mapping with spectra and lightcurves.

"""

import numpy as np
from LightCurve import LightCurve
from CCF import CCF

def generateLightCurve(spectra, minimum, maximum):
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

def generateLightCurveHBeta(spectra):
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

def getCrossCorrelationFunction(lightcurve1, lightcurve2, times=np.arange(-50., 100., 1.)):
    """
    Returns a cross correlation function (CCF) instance of a CCF class
    """

    # Generate lightcurves
    #self.generateLightCurve(5183, 5193, label='Cont')
    #self.generateLightCurveHBeta('wide', label='HBeta')

    # Get time coordinates used to generate the CCF
    lcCont = lightcurve1
    lcHBeta = lightcurve2

    # # # Do cross-correlation
    ccf1 = np.zeros(len(times))
    ccf2 = np.zeros(len(times))

    # Compute CCF1 where
    # each emission-line measurement L(t_i)
    # is paired with interpolated continuum values at C(t_i - tau)
    for i, tau in enumerate(times):
        for j in range(len(lcCont.time)):
            lagTime = lcCont.time[j] - tau
            if lagTime < min(lcCont.time) or lagTime > max(lcCont.time): continue
            ccf1[i] += (lcHBeta.flux[j] - np.mean(lcHBeta.flux)) \
                   * (lcCont.getFluxInterpolated(lagTime) - np.mean(lcCont.flux)) \
                   / (np.std(lcHBeta.flux) * np.std(lcCont.flux))
        ccf1[i] /= len(lcCont.time)

    # Compute CCF2 where
    # the measured continuum points C(t_i)
    # are paired with interpolated emission-line values at L(t_i + tau)
    for i, tau in enumerate(times):
        for j in range(len(lcCont.time)):
            lagTime = lcCont.time[j] + tau
            if lagTime < min(lcCont.time) or lagTime > max(lcCont.time): continue
            ccf2[i] += (lcCont.flux[j] - np.mean(lcCont.flux)) \
                   * (lcHBeta.getFluxInterpolated(lagTime) - np.mean(lcHBeta.flux)) \
                   / (np.std(lcHBeta.flux) * np.std(lcCont.flux))
        ccf2[i] /= len(lcCont.time)

    # Calculate mean CCF
    ccf = (ccf1 + ccf2) / 2.

    return CCF(times, ccf)

# def plotLightCurves(lightcurves):
#     """
#     Plot a range of lightcurves on top of eachother
#     """
#     # Generate lightcurves
#     self.generateLightCurve(5183, 5193, label='Cont')
#     self.generateLightCurveHBeta('wide', label='HBeta')
#
#     for i, lc in enumerate(self.lightcurves):
#         time = self.lightcurves[lc].time
#         flux = self.lightcurves[lc].flux
#         ferr = self.lightcurves[lc].ferr
#
#         plotID = int(str(len(self.lightcurves))+'1'+str(i))
#         plt.subplot(plotID)
#         plt.xlabel('time')
#         plt.ylabel('flux')
#         plt.legend(frameon=False)
#
#         plt.errorbar(time, flux, ferr, fmt='o', label=lc)
#
#     plt.title('Light curves from spectra')
#     plt.show()