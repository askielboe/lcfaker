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

    def getCrossCorrelationFunction(self, times=np.arange(-50., 100., 1.)):
        """
        Returns a cross correlation function (CCF) instance of a CCF class
        """

        # Generate lightcurves
        #self.generateLightCurve(5183, 5193, label='Cont')
        #self.generateLightCurveHBeta('wide', label='HBeta')

        # Get time coordinates used to generate the CCF
        lcCont = self.lcCont
        lcLine = self.lcLine

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
                ccf1[i] += (lcLine.flux[j] - np.mean(lcLine.flux)) \
                       * (lcCont.getFluxInterpolated(lagTime) - np.mean(lcCont.flux)) \
                       / (np.std(lcLine.flux) * np.std(lcCont.flux))
        ccf1 /= len(lcCont.time)

        # Compute CCF2 where
        # the measured continuum points C(t_i)
        # are paired with interpolated emission-line values at L(t_i + tau)
        for i, tau in enumerate(times):
            for j in range(len(lcCont.time)):
                lagTime = lcCont.time[j] + tau
                if lagTime < min(lcCont.time) or lagTime > max(lcCont.time): continue
                ccf2[i] += (lcCont.flux[j] - np.mean(lcCont.flux)) \
                       * (lcLine.getFluxInterpolated(lagTime) - np.mean(lcLine.flux)) \
                       / (np.std(lcLine.flux) * np.std(lcCont.flux))
        ccf2 /= len(lcCont.time)

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

        # Create mask of points in common betweent he lightcurves
        maskCont = self.lcCont.time >= min(self.lcLine.time)
        maskLine = self.lcLine.time <= max(self.lcCont.time)

        #
        self.lcCont.time = self.lcCont.time[maskCont]
        self.lcCont.flux = self.lcCont.flux[maskCont]

        self.lcLine.time = self.lcLine.time[maskLine]
        self.lcLine.flux = self.lcLine.flux[maskLine]

        if len(self.lcLine.ferr) > 0:
            self.lcCont.ferr = self.lcCont.ferr[maskCont]
            self.lcLine.ferr = self.lcLine.ferr[maskLine]

        # Renormalize time
        self.lcCont.time -= min(self.lcCont.time)
        self.lcLine.time -= min(self.lcLine.time)

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