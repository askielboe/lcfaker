class FakeLight():
	import numpy as np
	
	timeCont = np.array([])
	fluxCont = np.array([])
	timeLine = np.array([])
	fluxLine = np.array([])
	
	def __init__(self, length=1000):
		
		# Input parameters
		self.length = length
		
		# Generate a lightcurve
		from lcfaker import LightCurveMacLeod
		self.lightCurve = LightCurveMacLeod()
		
		from copy import copy
		self.lightCurveCont = copy(self.lightCurve)
		self.lightCurveLine = copy(self.lightCurve)
		
	def restore(self):
		from copy import copy
		self.lightCurveCont = copy(self.lightCurve)
		self.lightCurveLine = copy(self.lightCurve)
	
	def rebin(self, nBins):
		self.lightCurveCont.rebin(nBins)
		self.lightCurveLine.rebin(nBins)
	
	def reprocess(self, nBins, scale, sigma, c, alpha, beta):
		
		self.restore()
		
		# Reprocess the continuum to produce the line
		self.lightCurveLine.lag_luminosity(self.lightCurveCont, c, alpha, beta)
		
		import numpy as np
		import lib.units as units
		# print(max(self.lightCurveCont.flux))
		# print(units.mag_to_lum5100(max(self.lightCurveCont.flux)))
		# 
		# minLag = units.r_from_l(units.mag_to_lum5100(max(self.lightCurveCont.flux)))
		# maxLag = units.r_from_l(units.mag_to_lum5100(min(self.lightCurveCont.flux)))
		# avgLag = units.r_from_l(units.mag_to_lum5100(np.mean(self.lightCurveCont.flux)))
		# 
		# print "Minimum lag = ", minLag
		# print "Maximum lag = ", maxLag
		# print "Lag difference = ", maxLag - minLag
		# print "Lag corresponding to average luminosity (Continuum) = ", avgLag
		
		self.rebin(nBins)
		
		self.trim()
		
		self.lightCurveLine.scale(scale)
		self.lightCurveLine.smooth(sigma)
		
		return (minLag,maxLag,avgLag)
		
	def plot(self, style='-'):
		import matplotlib.pyplot as plt
		plt.figure()
		plt.plot(self.lightCurveCont.time,self.lightCurveCont.flux, str(style)+'b', label='cont')
		plt.plot(self.lightCurveLine.time,self.lightCurveLine.flux, str(style)+'r', label='line')
		plt.legend(frameon=False)
	
	def trim(self):
		# Find differences
		lowerDiff = int(min(self.lightCurveLine.time) - min(self.lightCurveCont.time))
		upperDiff = int(max(self.lightCurveLine.time) - max(self.lightCurveCont.time))
		
		# Remove first lowerDiff elements from lightCurveCont
		self.lightCurveCont.time = self.lightCurveCont.time[lowerDiff:]
		self.lightCurveCont.flux = self.lightCurveCont.flux[lowerDiff:]
		
		# Remove last upperDiff elements from lightCurveLine
		self.lightCurveLine.time = self.lightCurveLine.time[:len(self.lightCurveLine.time)-upperDiff]
		self.lightCurveLine.flux = self.lightCurveLine.flux[:len(self.lightCurveLine.flux)-upperDiff]
		
	def observeIntervals(self, times, widths):
		
		timeObserved = []
		fluxContObserved = []
		fluxLineObserved = []
		
		for i in range(len(times)):
			for j in range(len(self.lightCurveCont.time)):
				if self.lightCurveCont.time[j] > times[i]-widths[i] and self.lightCurveCont.time[j] < times[i]+widths[i]:
					timeObserved.append(self.lightCurveCont.time[j])
					fluxContObserved.append(self.lightCurveCont.flux[j])
					fluxLineObserved.append(self.lightCurveLine.flux[j])
		
		from numpy import asarray
		self.lightCurveCont.time = asarray(timeObserved)
		self.lightCurveLine.time = asarray(timeObserved)
		
		self.lightCurveCont.flux = asarray(fluxContObserved)
		self.lightCurveLine.flux = asarray(fluxLineObserved)
		
		#self.plot('o')
	
	def observeRandom(self, nRuns, lengthRun):
		
		# Generate observing runs randomly
		import random
		
		xsIntersection = list(set(self.lightCurveCont.time) & set(self.lightCurveLine.time))
		
		randomPositions = []
		for i in range(nRuns):
			randomPositions.append(random.choice(xsIntersection))
		
		xsObserved = []
		ysContObserved = []
		ysLineObserved = []
		
		xs = self.lightCurveCont.time
		ysCont = self.lightCurveCont.flux
		ysLine = self.lightCurveLine.flux
		
		# Pick out the given bins
		for position in randomPositions:
			for i in range(len(xs)):
				if xs[i] > position-lengthRun and xs[i] < position+lengthRun:
					xsObserved.append(xs[i])
					ysContObserved.append(ysCont[i])
					ysLineObserved.append(ysLine[i])
		
		self.lightCurveCont.time = xsObserved
		self.lightCurveLine.time = xsObserved
		
		self.lightCurveCont.flux = ysContObserved
		self.lightCurveLine.flux = ysLineObserved
		
		#self.plot('o')
	
	def observeCampaign(self, observingTime, numberOfObservations, scatterOfObservingTime):
		timeObserved = []
		fluxContObserved = []
		fluxLineObserved = []
		
		# observingTime: [0-1] should be related to AGN variability timescale
		# for now length is the total fractional observing time per lightcurve
		
		# numberOfObservations: [0-infinity] approximate number of observations
		
		integratedObservationTime = int(observingTime*len(self.lightCurveCont.time))
		print("integratedObservationTime",integratedObservationTime)
		
		meanObservingDuration = integratedObservationTime/numberOfObservations
		sigmaObservingDuration = scatterOfObservingTime * meanObservingDuration
		
		meanGapDuration = (len(self.lightCurveCont.time) - integratedObservationTime) / numberOfObservations
		sigmaGapDuration = scatterOfObservingTime * meanGapDuration
		
		# Loop through lightcurve adding observations and gaps randomly
		import random
		t = 0
		remainingObservingTime = integratedObservationTime
		for i in range(len(self.lightCurveCont.time)):
			
			# Get a random observing duration
			randomObservingDuration = int(random.gauss(meanObservingDuration,sigmaObservingDuration))
			print("randomObservingDuration",randomObservingDuration)
			
			# Break loop if we have no observing time left
			if t+randomObservingDuration > len(self.lightCurveCont.time): break
			
			#print("t+randomObservingDuration",t+randomObservingDuration)
			
			# Do the observation
			for j in range(randomObservingDuration):
				timeObserved.append(self.lightCurveCont.time[t+j])
				fluxContObserved.append(self.lightCurveCont.flux[t+j])
				fluxLineObserved.append(self.lightCurveLine.flux[t+j])
			
			t += randomObservingDuration
			remainingObservingTime -= randomObservingDuration
			
			if remainingObservingTime < 0: break
			
			
			# Get a random gap duration
			randomGapDuration = int(random.gauss(meanGapDuration,sigmaGapDuration))
			print("randomGapDuration",randomGapDuration)
			
			# Break loop if we have no observing time left
			if t+randomGapDuration > len(self.lightCurveCont.time): break
			
			t += randomGapDuration
		
		from numpy import asarray
		self.lightCurveCont.time = asarray(timeObserved)
		self.lightCurveLine.time = asarray(timeObserved)
		
		self.lightCurveCont.flux = asarray(fluxContObserved)
		self.lightCurveLine.flux = asarray(fluxLineObserved)
	
	def saveToTxt(self,fnamecont="lc_cont.txt",fnameline="lc_line.txt"):
		self.lightCurveCont.saveToTxt(fnamecont)
		self.lightCurveLine.saveToTxt(fnameline)
	
	def __repr__(self):
		print "<LightContAndLine(length='%f')" % (self.length)
		print "Average luminosity (Continuum) = '%f'" % self.lightCurveCont.getAverageFlux()
		print "Average luminosity (Line) = '%f'" % self.lightCurveLine.getAverageFlux()

