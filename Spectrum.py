# -*- coding: utf-8 -*-
"""
File:    Spectrum.py

Author:  Andreas Skielboe (skielboe@dark-cosmology.dk)
Date:    March 2013

Summary:

    Spectrum class.

"""

import numpy as np
import matplotlib.pyplot as plt

class Spectrum():
    """
    Class for loading, plotting and manipulating spectra
    Note: File names are n5xxxxyy.spc, where xxxx = four least significant
    figures in the Julian Date, and yy is a code indicating the origin of the data,
    given in the original papers.
    """

    def __init__(self):
        self.date = int(0)
        self.wavelength = np.array([])
        self.flux = np.array([])
        self.fluxErr = np.array([])

    def mask(self, minimum, maximum=-1):
        if maximum == -1:
            maximum = minimum
        mask = self.wavelength >= minimum
        mask *= self.wavelength <= maximum
        return mask

    def plot(self):
        plt.figure()

        # Plot errorbars if we have them
        if len(self.fluxErr) > 0:
            plt.errorbar(self.wavelength, self.flux, self.fluxErr)
        else:
            plt.plot(self.wavelength, self.flux)

        plt.title('spectrum')
        plt.xlabel('wavelength')
        plt.ylabel('flux')

        plt.show()

    # def loadData(self, filename):
    #     """
    #     Method for loading data
    #     (specialized for 'Selected Revised Optical Spectra of NGC 5548')
    #     http://www.astronomy.ohio-state.edu/~agnwatch/n5548/spectra/Revised/
    #     See: TODO.txt
    #     """
    #     self.name = filename.split('/')[-1]
    #     self.date = self.name[2:6]
    #
    #     data = np.loadtxt(filename)
    #
    #     self.wavelength = data.T[0]
    #     self.flux = data.T[1]
    #     self.fluxError = data.T[2]

    def integrate(self, minWavelength, maxWavelength):
        """
        Simple Riemann integration, summing bins.
        """
        mask = self.mask(minWavelength, maxWavelength)

        flux = self.flux[mask]
        fluxError = self.fluxError[mask]

        integralFlux = np.sum(flux)
        integralError = np.sqrt(np.sum(fluxError**2.))

        integralFlux *= 1./len(flux)
        integralError *= 1./len(flux)

        return integralFlux, integralError

    def fitContinuum(self, windows):
        """
        Fits the continuum using a 2nd degree polynomial using wavelength windows.
        Input: List of tuples defining the wavelength windows.
        Example: fitContinuumInWindows([(4600,4740),(4790,4840),(5130,5300)])
        """
        # Pick out the values we want to fit, from the complete spectrum
        mask = self.mask(-1)
        for window in windows:
            mask += self.mask(window[0],window[1])

        wavelength = self.wavelength[mask]
        flux = self.flux[mask]
        fluxError = self.fluxError[mask]

        # Fit a polynomial
        fit = np.polyfit(wavelength, flux, 2)

        # Calculate error from the continuum fit
        # using the average error in the bins used
        # divided by sqrt(numberOfBins)
        fitError = np.mean(fluxError) / np.sqrt(len(flux))

        return fit, fitError

    def integrateLine(self, window, continuum=-1, continuumErr=-1):
        """
        Subtracts continuum using straight line fit to two 20 A wide windows
        centered on 4840 and 5160.

        The total broad HB flux time series is obtained by summing
        all flux above an interpolated continuum in
        the observed wavelength band 4870-5010 A.

        TODO: Make more general line integration method,
        where the user defines integration limits, etc.
        """

        # Set limits for HBeta line
        mask = self.mask(4870, 5010)

        # Get parameter for the linear fit
        fit, fitError = self.fitContinuum([(4600,4740),(4790,4840),(5130,5300)])

        # Get polynomial based on fit parameters
        p = np.poly1d(fit)

        # Limit data to the HBeta line
        wavelength = self.wavelength[mask]
        flux = self.flux[mask]
        fluxError = self.fluxError[mask]

        # Calculate the integral flux with continuum subtracted
        integralFlux = np.sum(flux - p(wavelength))

        # Calculate the error for the line integration
        integralError = np.sqrt(np.sum(fluxError**2.) + len(flux)*(fitError**2.))

        # Normalize to window size (getting flux per 2 AA) - Unnecessary?
        integralFlux *= 1./len(flux)
        integralError *= 1./len(flux)

        return integralFlux, integralError

    def integrateHBeta(self):
        """
        Subtracts continuum using straight line fit to two 20 A wide windows
        centered on 4840 and 5160.

        The total broad HB flux time series is obtained by summing
        all flux above an interpolated continuum in
        the observed wavelength band 4870-5010 A.

        TODO: Make more general line integration method,
        where the user defines integration limits, etc.
        """

        # Set limits for HBeta line
        mask = self.mask(4870, 5010)

        # Get parameter for the linear fit
        fit, fitError = self.fitContinuum([(4600,4740),(4790,4840),(5130,5300)])

        # Get polynomial based on fit parameters
        p = np.poly1d(fit)

        # Limit data to the HBeta line
        wavelength = self.wavelength[mask]
        flux = self.flux[mask]
        fluxError = self.fluxError[mask]

        # Calculate the integral flux with continuum subtracted
        integralFlux = np.sum(flux - p(wavelength))

        # Calculate the error for the line integration
        integralError = np.sqrt(np.sum(fluxError**2.) + len(flux)*(fitError**2.))

        # Normalize to window size (getting flux per 2 AA) - Unnecessary?
        integralFlux *= 1./len(flux)
        integralError *= 1./len(flux)

        return integralFlux, integralError

    def addNoiseGaussian(self, snr):
        """
        Adding gaussian noise to the spectrum based on the provied SNR
        Note: The SNR is then per bin (e.g. 2 A)
        """
        squaredError = self.flux**2. / snr**2. - self.fluxError**2.
        if (squaredError < 0.).any():
            print """WARNING in addNoiseGaussian: Signal to noise larger than
                the instrinsic for the spectrum. No noise will be added!"""
            return

        std_noise = np.sqrt(squaredError)
        noise = std_noise * np.random.randn(len(self.flux))

        self.flux += noise
        self.fluxError = np.sqrt(self.fluxError**2. + std_noise**2.)

    def fit(self, minWavelength, maxWavelength, doplot=0):
        """
        LEGACY
        """
        x = self.wavelength[self.mask(minWavelength, maxWavelength)]
        y = self.flux[self.mask(minWavelength, maxWavelength)]

        def fitfunc(x, p):
            a, b, n, m1, m2, m3, s1, s2, s3 = p
            fitfunc = n*np.exp(-(x-m1)**2. / (2*s1**2.)) \
                + n*np.exp(-(x-m2)**2. / (2*s2**2.)) \
                + n*np.exp(-(x-m3)**2. / (2*s3**2.)) \
                + a*x + b
            return fitfunc

        def residuals(p, x, y):
            a, b, n, m1, m2, m3, s1, s2, s3 = p
            err = y - fitfunc(x, p)
            return err

        # Set starting parameters
        p0 = np.zeros(9)
        p0[0] = 0.     # a
        p0[1] = 0.     # b
        p0[2] = 0.76e-14   # n
        p0[3] = 4943.  # m1
        p0[4] = 4930.  # m2
        p0[5] = 4930.  # m3
        p0[6] = 7.     # s1
        p0[7] = 40.    # s2
        p0[8] = 50.    # s3

        from scipy.optimize import leastsq
        plsq = leastsq(residuals, p0, args=(x, y), full_output=0)

        def plot():
            plt.figure()
            plt.plot(x, fitfunc(x, plsq[0]), x, y, 'o')
            plt.title('Least-squares fit to spectrum')
            plt.legend(['Fit', 'Data'])
            plt.show()

        if doplot != 0:
            plot()

        return plsq

    def __repr__(self):
        return "Spectrum Instance: Spectrum wavelength range: " \
            + str(min(self.wavelength)) + ' - ' + str(max(self.wavelength))

