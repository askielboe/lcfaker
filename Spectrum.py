# -*- coding: utf-8 -*-
"""
File:    Spectrum.py

Author:  Andreas Skielboe (skielboe@dark-cosmology.dk)
Date:    March 2013

Summary:

    Spectrum class.

"""

from time import clock
import numpy as np
import matplotlib.pyplot as plt

class Spectrum():
    """
    Class for loading, plotting and manipulating spectra
    Note: File names are n5xxxxyy.spc, where xxxx = four least significant
    figures in the Julian Date, and yy is a code indicating the origin of the data,
    given in the original papers.
    """

    def __init__(self, wavelength = np.array([]), flux = np.array([]), ferr = np.array([]), date = int(-1)):
        self.date = date
        self.wavelength = wavelength
        self.flux = flux
        self.ferr = ferr

        # Parameters specific to emission line diagnostics
        # FIXME: Add support for multiple emission/absorption lines
        self.snrWindow = ()
        self.lineIntWindow = ()
        self.contIntWindow = ()
        self.contFitWindow = [()]
        self.oiiiWindow = ()

        if len(self.ferr) == 0:
            self.ferr = np.zeros(len(self.wavelength))

    def mask(self, minimum, maximum=-1):
        if maximum == -1:
            maximum = minimum
        mask = self.wavelength >= minimum
        mask *= self.wavelength <= maximum
        return mask

    def plot(self, showWindows=1, showContFit=1, saveFig=0):
        plt.figure()

        if showWindows:
            # Plot integration windows
            # FIXME: Use error handler instead of in-line exceptions!
            self.errorHandle("windows")

            plt.plot([self.lineIntWindow[0],self.lineIntWindow[0]], [min(self.flux),max(self.flux)], color='r')
            plt.plot([self.lineIntWindow[1],self.lineIntWindow[1]], [min(self.flux),max(self.flux)], color='r')

            for window in self.contFitWindow:
                plt.plot([window[0],window[0]], [min(self.flux),max(self.flux)], color='g')
                plt.plot([window[1],window[1]], [min(self.flux),max(self.flux)], color='g')

        if showContFit:
            self.errorHandle("windows")

            fitMean, fitCov = self.fitContinuum(self.contFitWindow)

            x = np.arange(min(self.wavelength), max(self.wavelength), 1.0)
            y = fitMean[0] * x + fitMean[1]

            plt.plot(x,y)

        # Plot errorbars if we have them
        if len(self.ferr) > 0:
            plt.errorbar(self.wavelength, self.flux, self.ferr)
        else:
            plt.plot(self.wavelength, self.flux)

        plt.title('spectrum')
        plt.xlabel('wavelength')
        plt.ylabel('flux')

        if saveFig:
            plt.xlim((4550,5450))
            plt.ylim((0.0, 1.0e-13))
            plt.savefig('tmp/spec_'+str(self.date)+'.jpg')
        else:
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

    #     data = np.loadtxt(filename)

    #     self.wavelength = data.T[0]
    #     self.flux = data.T[1]
    #     self.ferr = data.T[2]

    def getSNR(self, minWavelength = -1, maxWavelength = -1):
        """
        Return estimated signal to noise per spectral bin,
        assuming a flat continuum within given wavelength ranges
        Input: Wavelength range given by min, max of flat continuum
        Output: Estimated signal to noise per spectral bin
        """

        def getMaskedData(minWavelength, maxWavelength):
            mask = self.mask(minWavelength, maxWavelength)
            wavelength = self.wavelength[mask]
            flux = self.flux[mask]
            ferr = self.ferr[mask]
            return wavelength, flux, ferr

        if minWavelength == -1 or maxWavelength == -1:
            minWavelength = self.snrWindow[0]
            maxWavelength = self.snrWindow[1]

        wavelength, flux, ferr = getMaskedData(minWavelength, maxWavelength)

        # If the resolution is low, or the window too narrow, post a warning and increase window size
        while len(wavelength) < 3:
            print "\n Spectrum.getSNR: WARNING: Resolution / WindowSize too low! Increasing window size.. \n"
            minWavelength = minWavelength - (maxWavelength - minWavelength) / 2
            maxWavelength = maxWavelength + (maxWavelength - minWavelength) / 2
            wavelength, flux, ferr = getMaskedData(minWavelength, maxWavelength)

        # Remove possible slope from data by subtracting a straight line fit
        fit = np.polyfit(wavelength, flux, 1)
        p = np.poly1d(fit)

        slope = fit[0]

        # Rotate data about the pivot point, as to remove the slope
        # If we have an even number of data points the pivot point is between two datapoints
        # and we should take this into account
        if len(wavelength) % 2 == 0:
            iLeft = len(wavelength)/2
            iRight = len(wavelength)/2 + 1
            pivotFlux = p((wavelength[iLeft] + wavelength[iRight]) / 2.)

            flux[:iLeft] += np.sign(slope) * np.abs(p(wavelength[:iLeft]) - pivotFlux)
            flux[iRight:] -= np.sign(slope) * np.abs(p(wavelength[iRight:]) - pivotFlux)

        # If we have an uneven number of datapoints, the pivot point is on top of a datapoint,
        # and this point should not be altered
        elif len(wavelength) % 2 > 0:
            iLeft = len(wavelength)/2
            iRight = len(wavelength)/2 + 2
            iMiddle = len(wavelength)/2 + 1
            pivotFlux = p(wavelength[iMiddle])

            flux[:iLeft] += np.sign(slope) * np.abs(p(wavelength[:iLeft]) - pivotFlux)
            flux[iRight:] -= np.sign(slope) * np.abs(p(wavelength[iRight:]) - pivotFlux)

        snr = np.mean(flux)/np.std(flux)

        return snr

        # """
        # Simple Riemann integration, summing bins.
        # """
        # mask = self.mask(minWavelength, maxWavelength)

        # flux = self.flux[mask]
        # ferr = self.ferr[mask]

        # integralFlux = np.sum(flux)
        # integralError = np.sqrt(np.sum(ferr**2.))

        # integralFlux *= 1./len(flux)
        # integralError *= 1./len(flux)

        # return integralFlux, integralError

    def fitContinuum(self, contFitWindow=[()], plot=False):
        """
        Fits the continuum using a straight line using wavelength windows.
        Input: List of tuples defining the wavelength windows.
        Example: fitContinuumInWindows([(4600,4740),(4790,4840),(5130,5300)])
        Output: (fit parameters), (std of fit parameters based on cov matrix)
        """
        # Pick out the values we want to fit, from the complete spectrum
        mask = self.mask(-1)
        for window in contFitWindow:
            mask += self.mask(window[0],window[1])

        wavelength = self.wavelength[mask]
        flux = self.flux[mask]
        ferr = self.ferr[mask]

        # Fit a polynomial
        fitMean, fitCov = np.polyfit(wavelength, flux, 1, cov=True)

        if plot == True:
            plt.plot()

        return fitMean, fitCov

    def subtractContinuum(self, contFitWindow=[()]):
        fitMean, fitCov = self.fitContinuum(contFitWindow)
        continuum_fit = np.poly1d(fitMean)
        self.flux -= continuum_fit(self.wavelength)
        self.ferr *= np.sqrt(2)

    def integrateContinuum(self, contIntWindow=(), getErrors = 1, verbose = 0):
        """
        Calculates the flux integral of the given wavelength range
        without subtracting any continuum
        """
        if contIntWindow == ():
            contIntWindow = self.contIntWindow

        self.errorHandle("contIntWindow")

        # Pick out the values we want to fit, from the complete spectrum
        mask = self.mask(contIntWindow[0],contIntWindow[1])

        wavelength = self.wavelength[mask]
        flux = self.flux[mask]
        ferr = self.ferr[mask]

        aangstorms = np.max(wavelength) - np.min(wavelength)
        integralFluxCont = np.sum(flux) / aangstorms
        integralErrorCont = np.sqrt(np.sum(ferr**2.)) / aangstorms

        if verbose:
            print "continuum flux = ", integralFluxCont, " +/- ", integralErrorCont
            print "relative error cont flux (percent) = ", integralErrorCont/integralFluxCont * 100.0
            print "=================================================================="

        return integralFluxCont, integralErrorCont

    def integrateLine(self, getErrors = 1, verbose = 0, weight='uniform'):
        """
        Calculated integral flux by subtracting linearly interpolated continuum.
        Returns line flux and errer as well as estimated continuum flux and error.
        """

        self.errorHandle("windows")

        # Set limits for the line
        mask = self.mask(self.lineIntWindow[0],self.lineIntWindow[1])

        # Get parameter for the linear fit
        # EG: contFitWindow = [(4600,4740),(4790,4840),(5130,5300)]
        fitMean, fitCov = self.fitContinuum(self.contFitWindow)
        if verbose:
            print "cont fitMean = ", fitMean
            print "cont fitCov = ", fitCov

        # Get polynomial based on fit parameters
        p = np.poly1d(fitMean)

        # Limit data to the HBeta line
        wavelength = self.wavelength[mask]
        flux = self.flux[mask]
        ferr = self.ferr[mask]

        # Calculate the integral continuum
        # Remember the fit only makes sense at the wavelength bins given in the spectrum,
        # as non-constant bin widths are fitted as well.
        # Beware: THE FIT SHOULD ONLY BE EVALUATED AT BIN CENTERS GIVEN IN THE SPECTRUM
        # NO INTEGRALS, ETC.!
        integralFluxCont = np.sum(p(wavelength))

        # Calculate the integral line flux with continuum subtracted
        # Weigh each wavelength bin according to Horne (1986)' optimal extraction algorithm
        if weight == 'horne':
            # Weights based on Horne (1986)
            fractionalFlux = flux / np.sum(flux)
            w = (fractionalFlux / ferr)**2.
            flux = np.sum(fractionalFlux * flux / ferr**2.) / np.sum(w)
            ferr = np.sqrt(1./np.sum(w))
        elif weight == 'uniform':
            pass

        integralFluxLine = (np.sum(flux) - integralFluxCont)

        if getErrors == 1:
            # Calculate the error for the line integration using Monte Carlo
            # (avoids problems with uneven bin sizes)
            nSamples = 200

            # Draw fit parameters based on variance from the fit
            a, b = np.random.multivariate_normal(mean=fitMean, cov=fitCov, size=nSamples).T
            integrals = []
            for i in range(nSamples):
                integrals.append(np.sum(a[i] * wavelength + b[i]))

            if verbose:
                print "np.mean(integrals) = ", np.mean(integrals)
                print "np.std(integrals) = ", np.std(integrals)

            integralErrorCont = np.std(integrals)
            integralErrorLine = np.sqrt(np.sum(ferr**2.) + integralErrorCont**2.)
        else:
            integralErrorCont = 0.0
            integralErrorLine = 0.0

        if verbose:
            print "line flux = ", integralFluxLine, " +/- ", integralErrorLine
            print "continuum flux = ", integralFluxCont, " +/- ", integralErrorCont
            print "relative error line flux (percent) = ", integralErrorLine/integralFluxLine * 100.0
            print "relative error cont flux (percent) = ", integralErrorCont/integralFluxCont * 100.0
            print "=================================================================="

        return integralFluxLine, integralErrorLine, integralFluxCont, integralErrorCont

    def calibrateToOIII(self):
        # WORK IN PROGRESS!
        # Get polynomial based on fit parameters
        fitMean, fitCov = self.fitContinuum(self.contFitWindow)
        p = np.poly1d(fitMean)

        # Limit data to the OIII line
        mask = self.mask(self.oiiiWindow[0],self.oiiiWindow[1])
        wavelength = self.wavelength[mask]
        flux = self.flux[mask]
        ferr = self.ferr[mask]

        integralOIII = np.sum(flux) - np.sum(p(wavelength))

        self.flux /= integralOIII

    # TO BE DELETED!
    # def integrateHBeta(self):
    #     """
    #     NOTE: THIS ONLY WORKS ON WANDERS & PETERSON!
    #     Subtracts continuum using straight line fit to two 20 A wide windows
    #     centered on 4840 and 5160.

    #     The total broad HB flux time series is obtained by summing
    #     all flux above an interpolated continuum in
    #     the observed wavelength band 4870-5010 A.

    #     TODO: Make more general line integration method,
    #     where the user defines integration limits, etc.
    #     """

    #     integralFluxLine, integralErrorLine, integralFluxCont, integralErrorCont = self.integrateLine(
    #         lineIntWindow=(4870, 5010), contFitWindow=[(4600,4740),(4790,4840),(5130,5300)])

    #     # # Set limits for HBeta line
    #     # mask = self.mask(4870, 5010)

    #     # # Get parameter for the linear fit
    #     # fit, fitError = self.fitContinuum([(4600,4740),(4790,4840),(5130,5300)])

    #     # # Get polynomial based on fit parameters
    #     # p = np.poly1d(fit)

    #     # # Limit data to the HBeta line
    #     # wavelength = self.wavelength[mask]
    #     # flux = self.flux[mask]
    #     # ferr = self.ferr[mask]

    #     # # Calculate the integral flux with continuum subtracted
    #     # integralFlux = np.sum(flux - p(wavelength))

    #     # # Calculate the error for the line integration
    #     # integralError = np.sqrt(np.sum(ferr**2.) + len(flux)*(fitError**2.))

    #     # Normalize to window size (getting flux per 2 AA) - Unnecessary?
    #     #integralFlux *= 1./len(flux)
    #     #integralError *= 1./len(flux)

    #     return integralFluxLine, integralErrorLine

    def addNoiseGaussian(self, snr):
        """
        Adding gaussian noise to the spectrum based on the provied SNR
        Note: The SNR is then per bin (e.g. 2 A)
        Returns new Spectrum instance.
        """
        np.random.seed(int(clock()*10e6))

        squaredError = self.flux**2. / snr**2. - self.ferr**2.
        if (squaredError < 0.).any():
            print """WARNING in addNoiseGaussian: Signal to noise larger than
                the instrinsic for the spectrum. No noise will be added!"""
            return

        std_noise = np.sqrt(squaredError)
        noise = std_noise * np.random.randn(len(self.flux))

        #wavelength = self.wavelength.copy()
        self.flux = self.flux + noise
        self.ferr = np.sqrt(self.ferr**2. + std_noise**2.)
        #date = self.date

        #spectrum =


        #return Spectrum(wavelength = wavelength, flux = flux, ferr = ferr, date = date)

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

    def errorHandle(self, error):
        """
        Internal function for handling errors.
        """
        if error == "windows":
            if self.lineIntWindow == () or self.contIntWindow == () or self.contFitWindow == [()]:
                raise ValueError('Must define all integration windows (lineInt, contInt and contFit)!')
            if len(self.lineIntWindow) != 2:
                raise ValueError('lineIntWindow must have length == 2')

    def __repr__(self):
        return "Spectrum Instance: Spectrum wavelength range: " \
            + str(min(self.wavelength)) + ' - ' + str(max(self.wavelength))

