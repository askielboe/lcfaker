# -*- coding: utf-8 -*-
"""
File:    LightCurve.py

Author:  Andreas Skielboe (skielboe@dark-cosmology.dk)
Date:    March 2013

Summary:

    LightCurve class.

"""

import numpy as np

class LightCurve():

    def __init__(self, time=np.array([]), flux=np.array([]), ferr=np.array([]), label='lightcurve'):

        self.time = time
        self.flux = flux
        self.ferr = ferr
        self._label = label

    @property
    def label(self):
        return str(self._label)

    def mask(self, arr, minimum, maximum=-1):
        if maximum == -1:
            maximum = minimum
        mask = arr >= minimum
        mask *= arr <= maximum
        return mask

    def smooth(self, sigma=1.5):
        import numpy as np
        import scipy.signal as signal

        length = int(sigma*10)

        def gauss_kern(sigma=1.5):
            """ Returns a normalized 1D gauss kernel array for convolutions """
            x = np.array(range(length))
            x = x-length/2.
            g = np.exp( -( x**2. ) / (2.*sigma**2.) );
            return g / g.sum()

        g = gauss_kern(sigma)

        # To convolve close to boundaries we need to extrapolate beyond the boundaries
        # Constant extrapolation using only the values at the boundaries
        extraFluxLower = np.asarray([self.flux[0] for i in range(length)])
        extraFluxUpper = np.asarray([self.flux[-1] for i in range(length)])
        extendedFlux = np.concatenate((extraFluxLower,self.flux,extraFluxUpper))

        # Smooth light curve using numpy signal convolve and the gaussian kernel
        self.flux = np.asarray(signal.convolve(extendedFlux,g,mode='same'))
        self.flux = self.flux[length:]
        self.flux = self.flux[:-length]

    def lag_const(self, lag_const):
        print "Running lag_const.."

        self.time = self.time + lag_const

    # def lag_luminosity(self):
    #     """
    #     AGN Lightcurves: Method for lagging light curve based on a given continuum
    #     Input: Continuum lightcurve
    #     Effect: self will be modified to a lagged lightcurve.
    #     Output: None
    #
    #     TODO: Right now lag is based on R-L relation from Bentz et al.
    #     It only works if the light curve is generated from the MacLeod method/function.
    #     It _should_ work in general as a stand alone-method.
    #     """
    #     print "Running lag_luminosity.."
    #
    #     import numpy as np
    #     import lib.units as units
    #
    #     # Convert apparant magnitude to absolute magnitude before calculating the lag
    #     from copy import copy
    #     absoluteMag = copy(self.mag) - self.magApparant + self.magAbsolute
    #     print 'mean(absoluteMag)', np.mean(absoluteMag)
    #
    #     # Calculate timelag based on radius-luminosity relationship
    #     timelag = units.r_from_l(units.mag_to_lum5100(absoluteMag))
    #     self.time = self.time + timelag
    #
    #     # # #
    #     # Write some output
    #     print "Average luminosity: ", np.mean(units.mag_to_lum5100(absoluteMag))
    #
    #     minLag = units.r_from_l(units.mag_to_lum5100(max(absoluteMag)))
    #     maxLag = units.r_from_l(units.mag_to_lum5100(min(absoluteMag)))
    #     avgLag = units.r_from_l(units.mag_to_lum5100(np.mean(absoluteMag)))
    #
    #     print "Minimum lag = ", minLag
    #     print "Maximum lag = ", maxLag
    #     print "Lag difference = ", maxLag - minLag
    #     print "Lag corresponding to average luminosity (Continuum) = ", avgLag
    #     # # #
    #
    #     # Sort time array
    #     self.time.sort()

    def addNoiseGaussian(self, snr):
        """
        Add independent gaussian noise in each bin in the lightcurve.
        Usage: addNoiseGaussian([Signal To Noise])
        Output: None (modifies the instance).
        """
        if len(self.ferr) > 0:
            squaredError = self.flux**2. / snr**2. - self.ferr**2.
            if (squaredError < 0.).any():
                print """WARNING in addNoiseGaussian: Signal to noise larger than
                    the instrinsic for the lightcurve. No noise will be added!"""
                return

            std_noise = np.sqrt(squaredError)
            noise = std_noise * np.random.randn(len(self.flux))

            self.flux += noise
            self.ferr = np.sqrt(self.ferr**2. + std_noise**2.)

        elif len(self.ferr) == 0:
            std_noise = self.flux / snr
            noise = std_noise * np.random.randn(len(self.flux))

            self.flux += noise
            self.ferr = std_noise

    def rebin(self, shapeNew):
        """
        Rebin the lightcurve, summing fluxes and averaging times.
        Usage: rebin([new number of bins])
        Output: None (modifies the instance).
        Limitations: Only works if the current bin number is divisible by the new.
        """
        # First check is current array is divisible by nBins
        if (len(self.flux)%shapeNew > 0):
            print "ERROR: Current number of bins not divisible by new number of bins!"
            return

        # Rebin flux
        M = self.flux.shape[0]
        m = shapeNew
        self.flux = self.flux.reshape((m,M/m)).sum(1)

        # Calculate new time bins
        self.time = self.time.reshape((m,M/m)).mean(1)

    def scale(self, scale):
        """
        Rescales lightcurve by a constant (float).
        Output: None (modifies the instance).
        """
        self.flux = self.flux * float(scale)

    def getFluxInterpolated(self, t):
        """
        Returns a flux at any given time using
        linear interpolation between nearby data points

        t: time for which we calculate the interpolated flux
        f: interpolated flux at time t

        time: is an ndarray with J dates for the lightcurve
        flux: ndarray, the flux in the lightcurve
        """

        # First check if a datapoint exists at the given time
        mask = self.mask(self.time, t)
        if True in mask:
            f = self.flux[mask]
        else:
            maskLower = self.mask(self.time, 0, t)
            maskUpper = self.mask(self.time, t, max(self.time))

            indexLower = max(np.nonzero(maskLower)[0])
            indexUpper = min(np.nonzero(maskUpper)[0])

            # Calculate linear slope
            a = (self.flux[indexUpper] - self.flux[indexLower]) \
              / (self.time[indexUpper] - self.time[indexLower])

            # Calculated extrapolated flux at time t
            f = self.flux[indexLower] + a * (t - self.time[indexLower])

        return f

    def plot(self, units='Jy', figure=1, color='b', label='', marker='-'):
        import matplotlib.pyplot as plt
        plt.figure(figure)
        plt.title('LightCurve')
        plt.xlabel('time [days]')

        if (units == 'Jy'):
            import lib.units as units
            y = self.flux
            plt.ylabel('flux [Jy]')
        elif (units == 'Mag'):
            y = self.mag
            plt.ylabel('flux [m_i]')
        else:
            print "ERROR: Unknown unit"
            return

        plt.plot(self.time, y, marker+color, label=label)
        plt.legend(frameon=False)
        #plt.savefig('lightcurve.pdf')
        plt.show()

    def getAverageFlux(self):
        return sum(self.flux)/len(self.flux)

    def saveToTxt(self, outFileName='lightcurve.txt', signalToNoise=10.):
        from numpy import savetxt, transpose
        errorBars = 1./float(signalToNoise)*self.flux
        savetxt(outFileName, transpose((self.time, self.flux, errorBars)))

