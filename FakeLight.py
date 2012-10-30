class FakeLight():
	
	timeCont = []
	fluxCont = []
	timeLine = []
	fluxLine = []
	
	def __init__(self, length=1000):
		
		# Input parameters
		self.length = length
		# self.nBins = nBins
		# self.scale = scale
		# self.sigma = sigma
		# self.alpha = alpha
		# self.beta = beta
		
		# Generate a lightcurve
		import LightCurve
		self.lightCurve = LightCurve.LightCurve(self.length)
		
		# Now done when defining the light curve
		#lightCurve.generateWiener()
		#lightCurve.generateOU()
		
		# Normalize
		
		# for i in range(len(self.lightCurve.flux)):
			#self.lightCurve.flux[i] = self.lightCurve.flux[i] / max(self.lightCurve.flux) * 5. + 5.
		
		# Normalize using numpy array
		self.lightCurve.path[:,1] = (self.lightCurve.path[:,1] / max(abs(self.lightCurve.path[:,1])) + 1.)*5.
		
		from copy import copy
		self.lightCurveCont = copy(self.lightCurve)
		self.lightCurveLine = copy(self.lightCurve)
		
		# Extract (time,flux) values into two seperate list pairs
		# self.timeCont = list(lightCurve.time)
		# self.fluxCont = list(lightCurve.flux)
		# self.timeLine = list(lightCurve.time)
		# self.fluxLine = list(lightCurve.flux)
		
		# # Make new lightCurves for Cont, Line and original for storing the original lightcurve
		# self.lightCurveCont = LightCurve(self.length)
		# self.lightCurveContOrg = LightCurve(self.length)
		# self.lightCurveLine = LightCurve(self.length)
		# self.lightCurveLineOrg = LightCurve(self.length)
		# 
		# # Save working lightcurves
		# self.lightCurveCont.path = asarray(list(lightCurve.path))
		# self.lightCurveLine.path = asarray(list(lightCurve.path))
		# 
		# # Save originals for possible restore
		# self.lightCurveContOrg.time = list(self.lightCurveCont.time)
		# self.lightCurveContOrg.flux = list(self.lightCurveCont.flux)
		# self.lightCurveLineOrg.time = list(self.lightCurveLine.time)
		# self.lightCurveLineOrg.flux = list(self.lightCurveLine.flux)
		
	def restore(self):
		# self.lightCurveCont.time = list(self.lightCurveContOrg.time)
		# self.lightCurveCont.flux = list(self.lightCurveContOrg.flux)
		# 
		# self.lightCurveLine.time = list(self.lightCurveLineOrg.time)
		# self.lightCurveLine.flux = list(self.lightCurveLineOrg.flux)
		
		from copy import copy
		self.lightCurveCont = copy(self.lightCurve)
		self.lightCurveLine = copy(self.lightCurve)
	
	def rebin(self, nBins):
		self.lightCurveCont.rebin(nBins)
		self.lightCurveLine.rebin(nBins)
	
	def reprocess(self, nBins, scale, sigma, c, alpha, beta):
		
		self.restore()
		
		# Reprocess the continuum to produce the line
		if (beta == 0. and alpha == 0.): ## REDUNDANT!
			# Constant time-lag
			self.lightCurveLine.lag_const(c)
		else:
			# Luminosity dependent time-lag
			self.lightCurveLine.lag_luminosity(self.lightCurveCont, c, alpha, beta)
			#self.lightCurveCont.bin(nBins)
			print "Lag corresponding to average luminosity (Continuum) = ", c + alpha * (self.lightCurveCont.getAverageFlux())**beta
			print "Lag corresponding to average luminosity (Line) = ", c + alpha * (self.lightCurveLine.getAverageFlux())**beta
			print "Minimum lag = ", c + alpha * (min(self.lightCurveCont.flux))**beta
			print "Maximum lag = ", c + alpha * (max(self.lightCurveCont.flux))**beta
			print "Lag difference = ", (c + alpha * (max(self.lightCurveCont.flux))**beta) - (c + alpha * (min(self.lightCurveCont.flux))**beta)
		
		self.trim()
		
		self.rebin(nBins)
		
		self.lightCurveLine.scale(scale)
		self.lightCurveLine.smooth(sigma)
		
		self.plot()
		
	def plot(self, style='-'):
		import matplotlib.pyplot as plt
		plt.figure()
		plt.plot(self.lightCurveCont.time,self.lightCurveCont.flux, str(style)+'b', label='cont')
		plt.plot(self.lightCurveLine.time,self.lightCurveLine.flux, str(style)+'r', label='line')
		plt.legend(frameon=False)
		#plt.show()
	
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
		
		self.plot('o')
	
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
		
		self.plot('o')
		
	def saveToTxt(self):
		self.lightCurveCont.saveToTxt('lc_cont.txt')
		self.lightCurveLine.saveToTxt('lc_line.txt')
	
	def __repr__(self):
		print "<LightContAndLine(length='%f')" % (self.length)
		print "Average luminosity (Continuum) = '%f'" % self.lightCurveCont.getAverageFlux()
		print "Average luminosity (Line) = '%f'" % self.lightCurveLine.getAverageFlux()

