# -*- coding: utf-8 -*-
"""
File:    lcfunctions.py

Author:  Andreas Skielboe (skielboe@dark-cosmology.dk)
Date:    March 2013

Summary:

    Functions for manipulating LightCurve class instances.

"""
import numpy as np
from LightCurve import LightCurve

def plot(lc1, lc2):
    """
    Function to plot two lightcurves
    """
    import matplotlib.pyplot as plt
    plt.figure()
    plt.errorbar(lc1.time, lc1.flux, lc1.ferr, label=lc1.label)
    plt.errorbar(lc2.time, lc2.flux, lc2.ferr, label=lc2.label)
    plt.show()

def trim(lc1, lc2):
    """
    Crops the lightcurves to only include the overlapping regions.
    Input: Two LightCurve instances
    Output: Two LightCurve instances including only overlapping regions of the input
    """

    lcout1 = LightCurve()
    lcout2 = LightCurve()

    # Create mask of points in common betweent he lightcurves
    maskCont = lc1.time >= min(lc2.time)
    maskLine = lc2.time <= max(lc1.time)

    #
    lcout1.time = lc1.time[maskCont]
    lcout1.flux = lc1.flux[maskCont]

    lcout2.time = lc2.time[maskLine]
    lcout2.flux = lc2.flux[maskLine]

    # Renormalize time
    lcout1.time -= min(lc1.time)
    lcout2.time -= min(lc2.time)

    return lcout1, lcout2

def observeConstantCadence(self, lc1, lc2, nObservations):
    """
    Function to 'observe' lightcurves with constant cadence.
    """

    lcout1 = LightCurve()
    lcout2 = LightCurve()

    timeObserved = []
    fluxContObserved = []
    fluxLineObserved = []

    fluxContError = []
    fluxLineError = []

    # Dividing integers will floor automagically
    cadence = int(len(lc1.time))/int(nObservations)

    for i in range(nObservations):
        t = cadence*i
        timeObserved.append(lc1.time[t])

        fluxContObserved.append(lc1.flux[t])
        fluxLineObserved.append(lc2.flux[t])

        # Adding 5 % errors
        fluxContError.append(0.05*lc1.flux[t])
        fluxLineError.append(0.05*lc2.flux[t])

    lcout1.time = np.asarray(timeObserved)
    lcout2.time = np.asarray(timeObserved)

    lcout1.flux = np.asarray(fluxContObserved)
    lcout2.flux = np.asarray(fluxLineObserved)

    lcout1.ferr = np.asarray(fluxContError)
    lcout2.ferr = np.asarray(fluxLineError)

    return lcout1, lcout2

# def observeIntervals(self, times, widths):
#
#     timeObserved = []
#     fluxContObserved = []
#     fluxLineObserved = []
#
#     for i in range(len(times)):
#         for j in range(len(self.lightCurveCont.time)):
#             if self.lightCurveCont.time[j] > times[i]-widths[i] and self.lightCurveCont.time[j] < times[i]+widths[i]:
#                 timeObserved.append(self.lightCurveCont.time[j])
#                 fluxContObserved.append(self.lightCurveCont.flux[j])
#                 fluxLineObserved.append(self.lightCurveLine.flux[j])
#
#     from numpy import asarray
#     self.lightCurveCont.time = asarray(timeObserved)
#     self.lightCurveLine.time = asarray(timeObserved)
#
#     self.lightCurveCont.flux = asarray(fluxContObserved)
#     self.lightCurveLine.flux = asarray(fluxLineObserved)
#
#     #self.plot('o')
#
# def observeRandom(self, nRuns, lengthRun):
#
#     # Generate observing runs randomly
#     import random
#
#     xsIntersection = list(set(self.lightCurveCont.time) & set(self.lightCurveLine.time))
#
#     randomPositions = []
#     for i in range(nRuns):
#         randomPositions.append(random.choice(xsIntersection))
#
#     xsObserved = []
#     ysContObserved = []
#     ysLineObserved = []
#
#     xs = self.lightCurveCont.time
#     ysCont = self.lightCurveCont.flux
#     ysLine = self.lightCurveLine.flux
#
#     # Pick out the given bins
#     for position in randomPositions:
#         for i in range(len(xs)):
#             if xs[i] > position-lengthRun and xs[i] < position+lengthRun:
#                 xsObserved.append(xs[i])
#                 ysContObserved.append(ysCont[i])
#                 ysLineObserved.append(ysLine[i])
#
#     self.lightCurveCont.time = xsObserved
#     self.lightCurveLine.time = xsObserved
#
#     self.lightCurveCont.flux = ysContObserved
#     self.lightCurveLine.flux = ysLineObserved
#
#     #self.plot('o')
#
# def observeCampaign(self, observingTime, numberOfObservations, scatterOfObservingTime):
#     timeObserved = []
#     fluxContObserved = []
#     fluxLineObserved = []
#
#     # observingTime: [0-1] should be related to AGN variability timescale
#     # for now length is the total fractional observing time per lightcurve
#
#     # numberOfObservations: [0-infinity] approximate number of observations
#
#     integratedObservationTime = int(observingTime*len(self.lightCurveCont.time))
#     print("integratedObservationTime",integratedObservationTime)
#
#     meanObservingDuration = integratedObservationTime/numberOfObservations
#     sigmaObservingDuration = scatterOfObservingTime * meanObservingDuration
#
#     meanGapDuration = (len(self.lightCurveCont.time) - integratedObservationTime) / numberOfObservations
#     sigmaGapDuration = scatterOfObservingTime * meanGapDuration
#
#     # Loop through lightcurve adding observations and gaps randomly
#     import random
#     t = 0
#     remainingObservingTime = integratedObservationTime
#     for i in range(len(self.lightCurveCont.time)):
#
#         # Get a random observing duration
#         randomObservingDuration = int(random.gauss(meanObservingDuration,sigmaObservingDuration))
#         print("randomObservingDuration",randomObservingDuration)
#
#         # Break loop if we have no observing time left
#         if t+randomObservingDuration > len(self.lightCurveCont.time): break
#
#         #print("t+randomObservingDuration",t+randomObservingDuration)
#
#         # Do the observation
#         for j in range(randomObservingDuration):
#             timeObserved.append(self.lightCurveCont.time[t+j])
#             fluxContObserved.append(self.lightCurveCont.flux[t+j])
#             fluxLineObserved.append(self.lightCurveLine.flux[t+j])
#
#         t += randomObservingDuration
#         remainingObservingTime -= randomObservingDuration
#
#         if remainingObservingTime < 0: break
#
#
#         # Get a random gap duration
#         randomGapDuration = int(random.gauss(meanGapDuration,sigmaGapDuration))
#         print("randomGapDuration",randomGapDuration)
#
#         # Break loop if we have no observing time left
#         if t+randomGapDuration > len(self.lightCurveCont.time): break
#
#         t += randomGapDuration
#
#     from numpy import asarray
#     self.lightCurveCont.time = asarray(timeObserved)
#     self.lightCurveLine.time = asarray(timeObserved)
#
#     self.lightCurveCont.flux = asarray(fluxContObserved)
#     self.lightCurveLine.flux = asarray(fluxLineObserved)
