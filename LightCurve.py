# -*- coding: utf-8 -*-
"""
File:    LightCurve.py

Author:  Andreas Skielboe (skielboe@dark-cosmology.dk)
Date:    March 2013

Summary:

    LightCurve class.

"""

import numpy as np
import matplotlib.pyplot as plt

class LightCurve():

    def __init__(self, time=np.array([]), flux=np.array([]), ferr=np.array([]), label='lightcurve'):

        # Sort values according to time
        self.flux = np.array([f for (t,f) in sorted(zip(time,flux))])
        self.ferr = np.array([fe for (t,fe) in sorted(zip(time,ferr))])
        self.time = np.sort(time)
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

    def smooth(self, sigma=3.0, mu = 0.0, type='gaus'):
        import numpy as np
        import scipy.signal as signal

        length = int(sigma*10)

        def gauss_kern(sigma):
            """ Returns a normalized 1D gauss kernel array for convolutions """
            x = np.array(range(length))
            x = x-length/2.
            g = np.exp( -( x**2. ) / (2.*sigma**2.) );
            return g / g.sum()

        def gauss_kern_truncated(mu, sigma):
            """ Returns a normalized 1D gauss kernel array for convolutions """
            x = np.arange(0., 25.1, 1.0)
            zeros = np.zeros(25)
            g = np.exp( -(mu - x)**2. / (2.*sigma**2.) );
            g = np.concatenate((zeros,g))

            return g / g.sum(), length

        def lorentz_kern(sigma):
            """ Returns a normalized 1D lorentz kernel array for convolutions """
            x = np.array(range(length))
            x = x-length/2.
            g = 1/np.pi * sigma / (x**2. + sigma**2.);
            return g / g.sum()

        def mrk50_transfer_function_1():
            """ Returns a normalized Mrk50 transfer function kernel """
            #norm = 31.3044
            a = 4.84766
            b = -0.392628
            c = -0.129766

            # Forcing f(0) = 0
            x = np.arange(1.0, 25.1, 1.0)
            y = a * np.exp(b * x ** -2.0 + c*x)
            y = np.concatenate((np.zeros(1),y))

            # Add zeros before t = 0 (causality)
            zeros = np.zeros(25)
            x = np.arange(-25.0, 25.1, 1.0)
            y = np.concatenate((zeros,y))

            length = len(x)
            return y/y.sum(), length

        def mrk50_transfer_function_2():
            """ Returns a normalized Mrk50 transfer function kernel """
            #norm = 33.3227
            a = 1.75
            b = 0.825
            c = 0.219

            x = np.arange(0., 25.1, 1.0)
            y = a * x**2. / (b + x*np.exp(c*x))

            # Add zeros before t = 0 (causality)
            zeros = np.zeros(25)
            x = np.arange(-25.0, 25.1, 1.0)
            y = np.concatenate((zeros,y))

            length = len(x)
            return y/y.sum(), length

        if type == 'gaus':
            g = gauss_kern(sigma)
        elif type == 'gaus_trunc':
            g, length = gauss_kern_truncated(mu, sigma)
        elif type == 'lorentz':
            g = lorentz_kern(sigma)
        elif type == 'mrk50_1':
            g, length = mrk50_transfer_function_1()
        elif type == 'mrk50_2':
            g, length = mrk50_transfer_function_2()

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
        #print "Running lag_const.."
        self.time = self.time + lag_const
        return None

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

        # If no error information is available, create new noise
        elif len(self.ferr) == 0 or True in (self.ferr != np.zeros(len(self.ferr))):
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

    def getFluxInterpolated(self, tInter):
        """
        Returns a flux at any given time using
        linear interpolation between nearby data points

        t: time for which we calculate the interpolated flux
        f: interpolated flux at time t

        time: is an ndarray with J dates for the lightcurve
        flux: ndarray, the flux in the lightcurve
        """

        # Make sure that tInter is a numpy array
        if type(tInter) != np.ndarray and type(tInter) != list:
            tInter = np.array([tInter])
        elif type(tInter) == list:
            tInter = np.array(tInter)

        # If chosen time is outside self time raise an exception
        if np.any(tInter < min(self.time)) or np.any(tInter > max(self.time)):
            raise ValueError('Time for interpolation outside data range. \n mintInt = '+str(min(tInter))+' maxtInt = '+str(max(tInter))+', \n minTime = '+str(min(self.time))+', maxTime = '+str(max(self.time)))

        # Otherwise we do interpolation
        # First calculate slope between all points
        # We do this by subtracting self.time array from a shifted verion of itself
        slopes = (self.flux[1:] - self.flux[:-1]) \
               / (self.time[1:] - self.time[:-1])

        # Calculate the fluxes we don't already have from data
        # Figure out which slopes belong to which tInter
        # BUT: If max(self.time) is in iInter we have a problem, so we take it out..
        maxIntInter = False
        if max(tInter) == max(self.time):
            tInter = tInter[:-1]
            maxIntInter = True

        index = []
        for ti in tInter:
            index.append(np.max(np.nonzero(self.time <= ti)))

        # The correct slopes for each tInter is then slopes[index]
        # And the interpolated fluxes are
        fInter = (tInter - self.time[index])*slopes[index] + self.flux[index]

        # If we have removed the last point from tInter, add this to the flux now
        if maxIntInter == True:
            fInter = np.append(fInter,self.flux[-1])

        return fInter

    def observeConstantCadence(self, nObs, snr=-1):
        """
        Function to 'observe' lightcurve with constant cadence.
        """

        if nObs > len(self.time):
            raise ValueError('Number of observations larger than available data! Please choose a lower nObs. \n \
                nObs = '+str(nObs)+' nDays = '+str(len(self.time)))

        lcObserved = LightCurve()

        timeObserved = []
        fluxObserved = []
        ferrObserved = []

        # Calculate observing times
        try:
            cadence = len(self.time)/float(nObs)
        except ZeroDivisionError:
            raise ZeroDivisionError("len(self.time) / float(nObs) = "+str(len(self.time))+"/"+str(float(nObs)))
        obsTimes = np.arange(0.,len(self.time)-1,cadence)

        # Add random shift to obsTimes
        #obsTimes += np.random.rand()*cadence
        #obsTimes = np.round(obsTimes)
        #obsTimes = obsTimes[obsTimes < len(self.time)]

        for t in obsTimes:
            timeObserved.append(self.time[t])
            fluxObserved.append(self.flux[t])
            ferrObserved.append(self.ferr[t])

        lcObserved.time = np.asarray(timeObserved)
        lcObserved.flux = np.asarray(fluxObserved)
        lcObserved.ferr = np.asarray(ferrObserved)

        # Adding noise
        if snr > -1:
            lcObserved.addNoiseGaussian(snr)

        return lcObserved

    def plotfft(self):
        fft = np.fft.fft(self.flux)
        ps = np.abs(fft)**2.

        time_step = 1. / 30.
        freqs = np.fft.fftfreq(self.flux.size, time_step)
        idx = np.argsort(freqs)

        print "max(ps) = ", max(ps)
        print "std(ps) = ", np.std(ps)
        plt.figure()
        plt.semilogy()
        plt.xlim(0,len(self.flux)/2)
        #plt.plot(p)
        plt.plot(freqs[idx], ps[idx])
        plt.show()

    def plot(self, units='Jy', figure=1, color='b', label='', marker='-'):
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

        if len(self.ferr) == len(self.flux):
            plt.errorbar(self.time, y, self.ferr, fmt=marker+color, label=label)
        else:
            plt.plot(self.time, y, marker=marker, label=label)

        plt.legend(frameon=False)
        #plt.savefig('lightcurve.pdf')
        plt.show()

    def getAverageFlux(self):
        return sum(self.flux)/len(self.flux)

    def saveToTxt(self, outFileName='lightcurve.txt', signalToNoise=10.):
        from numpy import savetxt, transpose
        errorBars = 1./float(signalToNoise)*self.flux
        savetxt(outFileName, transpose((self.time, self.flux, errorBars)))

