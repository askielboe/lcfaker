# -*- coding: utf-8 -*-
"""
File:    reverberation.py

Author:  Andreas Skielboe (skielboe@dark-cosmology.dk)
Date:    March 2013

Summary:

    Class for doing reverberation mapping using light curves.

"""

from copy import deepcopy
import numpy as np
from CCF import CCF
from LightCurve import LightCurve

class Reverberation():

    def __init__(self, lcCont, lcLine):
        # Deepcopy light curves instances, so they are not modified
        # Will-it-work!? - Memory leak!?
        # self.lcCont = deepcopy(lcCont)
        # self.lcLine = deepcopy(lcLine)

        self.lcCont = lcCont
        self.lcLine = lcLine

    def getCrossCorrelationFunction(self, minLag=-50, maxLag=100, resLag=1):
        """
        Returns a cross correlation function (CCF) instance of a CCF class.
        """

        # Get array of time shifts
        times = np.arange(minLag, maxLag, resLag)

        # Get light curves
        lcCont = self.lcCont
        lcLine = self.lcLine

        # Allocate memory
        ccf1 = np.zeros(len(times))
        ccf2 = np.zeros(len(times))

        # Get interpolated time and flux over the whole light curve range with resLag intervals
        # These are the fluxed used for the t_i +/- tau function calls in the CCF
        lcContTimeInterpolated = np.arange(min(lcCont.time), max(lcCont.time) + resLag, resLag)
        lcLineTimeInterpolated = np.arange(min(lcLine.time), max(lcLine.time) + resLag, resLag)
        lcContFluxInterpolated = lcCont.getFluxInterpolated(lcContTimeInterpolated)
        lcLineFluxInterpolated = lcLine.getFluxInterpolated(lcLineTimeInterpolated)

        # Compute CCF1 where
        # each emission-line measurement L(t_i)
        # is paired with interpolated continuum values at C(t_i - tau)
        for i, tau in enumerate(times):
            # Calculate the times for the lagged light curve
            timeLagged = lcLine.time - tau

            # Some of the lagged times we can't calculate, becuase they are outside data range
            # These we mask out in timeLagged
            timeLaggedMask = (timeLagged >= min(lcCont.time)) * (timeLagged <= max(lcCont.time))
            timeLagged = timeLagged[timeLaggedMask]
            if len(timeLagged) <= 1: continue

            # Now we need to find the fluxes for C(t_i - tau)
            # WARNING: This only works because the times are in increasing order!
            lcContTimeIndices = np.searchsorted(lcContTimeInterpolated, timeLagged)
            #lcLineTime = lcLine.time[timeLaggedMask] # We don't actually need these times as we can just put the mask in lcLine.flux directly (as done below). And the flux is all we need.
            # Here we set the actual fluxes
            lcContFluxInterLagged = lcContFluxInterpolated[lcContTimeIndices]
            # Since timeLagged is based on lcLine times, for each time - tau we can't calculate, the corresponding time
            # we remove in lcLine.
            lcLineFlux = lcLine.flux[timeLaggedMask]

            # Compute the CC value for the given tau
            ccf1[i] = np.sum((lcLineFlux - np.mean(lcLineFlux)) \
                   * (lcContFluxInterLagged - np.mean(lcContFluxInterLagged)) \
                   / (np.std(lcLineFlux) * np.std(lcContFluxInterLagged))) \
                   / len(lcLineFlux)

        # Compute CCF2 where
        # the measured continuum points C(t_i)
        # are paired with interpolated emission-line values at L(t_i + tau)
        for i, tau in enumerate(times):
            timeLagged = lcCont.time + tau

            # Some of the lagged times we can't calculate, becuase they are outside data range
            # These we mask out in timeLagged
            timeLaggedMask = (timeLagged >= min(lcLine.time)) * (timeLagged <= max(lcLine.time))
            timeLagged = timeLagged[timeLaggedMask]
            if len(timeLagged) <= 1: continue

            # Now we need to find the fluxes for L(t_i - tau)
            # WARNING: This only works because the times are in increasing order!
            lcLineTimeIndices = np.searchsorted(lcLineTimeInterpolated, timeLagged)

            # Here we set the actual fluxes
            lcLineFluxInterLagged = lcLineFluxInterpolated[lcLineTimeIndices]
            lcContFlux = lcCont.flux[timeLaggedMask]

            ccf2[i] = np.sum((lcContFlux - np.mean(lcContFlux)) \
                   * (lcLineFluxInterLagged - np.mean(lcLineFluxInterLagged)) \
                   / (np.std(lcContFlux) * np.std(lcLineFluxInterLagged))) \
                   / len(lcContFlux)

        # Calculate mean CCF
        ccf = (ccf1 + ccf2) / 2.

        return CCF(times, ccf)

    def plot(self, marker='o', linestyle=''):
        """
        Function to plot two lightcurves
        """
        import matplotlib.pyplot as plt
        plt.figure()
        if len(self.lcCont.ferr) > 0 and len(self.lcLine.ferr) > 0:
            plt.errorbar(self.lcCont.time, self.lcCont.flux, self.lcCont.ferr, label=self.lcCont.label, marker=marker, linestyle=linestyle)
            plt.errorbar(self.lcLine.time, self.lcLine.flux, self.lcLine.ferr, label=self.lcLine.label, marker=marker, linestyle=linestyle)
        else:
            plt.plot(self.lcCont.time, self.lcCont.flux, label=self.lcCont.label, marker=marker, linestyle=linestyle)
            plt.plot(self.lcLine.time, self.lcLine.flux, label=self.lcLine.label, marker=marker, linestyle=linestyle)
        plt.show()

    def trim(self):
        """
        Crops the lightcurves to only include the overlapping regions.
        Input: Two LightCurve instances
        Output: None (works on the instance).
        """

        # Create mask of points in common betweent the lightcurves
        maskCont = (self.lcCont.time >= min(self.lcLine.time)) * (self.lcCont.time <= max(self.lcLine.time))
        maskLine = (self.lcLine.time >= min(self.lcCont.time)) * (self.lcLine.time <= max(self.lcCont.time))

        self.lcCont.time = self.lcCont.time[maskCont]
        self.lcCont.flux = self.lcCont.flux[maskCont]

        self.lcLine.time = self.lcLine.time[maskLine]
        self.lcLine.flux = self.lcLine.flux[maskLine]

        if len(self.lcLine.ferr) > 0:
            self.lcCont.ferr = self.lcCont.ferr[maskCont]
            self.lcLine.ferr = self.lcLine.ferr[maskLine]

        # Renormalize time
        #self.lcCont.time -= min(self.lcCont.time)
        #self.lcLine.time -= min(self.lcLine.time)

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
