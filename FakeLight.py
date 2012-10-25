class LightContAndLine():
	
	timeCont = []
	fluxCont = []
	timeLine = []
	fluxLine = []
	
	def __init__(self, length, scale, sigma, alpha, beta):
		
		# Input parameters
		self.length = length
		self.scale = scale
		self.sigma = sigma
		self.alpha = alpha
		self.beta = beta
		
		# Generate a lightcurve
		lightCurve = LightCurveWiener(self.length)
		lightCurve.generate()
		
		# Extract (time,flux) values into two seperate list pairs
		# self.timeCont = list(lightCurve.time)
		# self.fluxCont = list(lightCurve.flux)
		# self.timeLine = list(lightCurve.time)
		# self.fluxLine = list(lightCurve.flux)
		
		# Make new lightCurves for Cont, Line and original for storing the original lightcurve
		self.lightCurveCont = LightCurveWiener(self.length)
		self.lightCurveContOrg = LightCurveWiener(self.length)
		self.lightCurveLine = LightCurveWiener(self.length)
		self.lightCurveLineOrg = LightCurveWiener(self.length)
		
		# Using the (time,flux) list pairs as extracted earlier
		self.lightCurveCont.time = list(lightCurve.time)
		self.lightCurveCont.flux = list(lightCurve.flux)
		
		self.lightCurveLine.time = list(lightCurve.time)
		self.lightCurveLine.flux = list(lightCurve.flux)
		
		# Reprocess the continuum
		self.lightCurveLine.smooth(self.sigma)
		self.lightCurveLine.scale(self.scale)
		
		# Constant time-lag
		self.lightCurveLine.lag_const(self.alpha)
		
		# Luminosity dependent time-lag
		# nBins = 100
		# self.lightCurveLine.lag_luminosity(nBins, self.alpha, self.beta)
		# self.lightCurveCont.bin(nBins)
		
		# Save originals for possible restore
		self.lightCurveContOrg.time = list(self.lightCurveCont.time)
		self.lightCurveContOrg.flux = list(self.lightCurveCont.flux)
		self.lightCurveLineOrg.time = list(self.lightCurveLine.time)
		self.lightCurveLineOrg.flux = list(self.lightCurveLine.flux)
	
	def restore(self):
		self.lightCurveCont.time = list(self.lightCurveContOrg.time)
		self.lightCurveCont.flux = list(self.lightCurveContOrg.flux)
		
		self.lightCurveLine.time = list(self.lightCurveLineOrg.time)
		self.lightCurveLine.flux = list(self.lightCurveLineOrg.flux)
	
	def bin(self, nBins):
		self.lightCurveCont.bin(nBins)
		self.lightCurveLine.bin(nBins)
	
	def plot(self, style='-'):
		import matplotlib.pyplot as plt
		plt.figure()
		plt.plot(self.lightCurveCont.time,self.lightCurveCont.flux, str(style)+'b', label='cont')
		plt.plot(self.lightCurveLine.time,self.lightCurveLine.flux, str(style)+'r', label='line')
		plt.legend(frameon=False)
		#plt.show()
	
	def doIntersection(self):
		#### BROKEN ####
		################
		# Cut away regions with no overlap
		xsIntersection = list(set(self.lightCurveCont.time) & set(self.lightCurveLine.time))
		
		ysContIntersection = []
		ysLineIntersection = []
		for x in xsIntersection:
			ysContIntersection.append(self.lightCurveCont.flux[int(x)])
			ysLineIntersection.append(self.lightCurveLine.flux[int(x)])
		
		self.lightCurveCont.time = xsIntersection
		self.lightCurveLine.time = xsIntersection
		
		self.lightCurveCont.flux = ysContIntersection
		self.lightCurveLine.flux = ysLineIntersection
	
	def trim(self):
		# Find differences
		lowerDiff = int(min(self.lightCurveLine.time) - min(self.lightCurveCont.time))
		upperDiff = int(max(self.lightCurveLine.time) - max(self.lightCurveCont.time))
		
		# Remove first lowerDiff elements from lightCurveCont
		self.lightCurveCont.time = list(self.lightCurveCont.time[lowerDiff:])
		self.lightCurveCont.flux = list(self.lightCurveCont.flux[lowerDiff:])
		
		# Remove last upperDiff elements from lightCurveLine
		self.lightCurveLine.time = list(self.lightCurveLine.time[:len(self.lightCurveLine.time)-upperDiff])
		self.lightCurveLine.flux = list(self.lightCurveLine.flux[:len(self.lightCurveLine.flux)-upperDiff])
		
	#def observeIntervals(self, positions, widths):
		
	
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
		
	def saveToTxt(self):
		self.lightCurveCont.saveToTxt('lc_cont.txt')
		self.lightCurveLine.saveToTxt('lc_line.txt')

class LightCurveWiener():
	
	time = []
	flux = []
	
	def __init__(self, length):
		
		self.length = length
		
		# Generate random walk light curve using pyprocesses.py
		self.wienerParams = {"mu":0, "sigma":1}
		self.wienerInitial = {"startTime":0, "startPosition":0, "endTime":1200, "endPosition": 1 }
		
	
	def generate(self):
		import pyprocess as SP
		Wiener = SP.Wiener_process(self.wienerParams, self.wienerInitial)
		path = Wiener.generate_sample_path(range(self.length))
		self.time = [x for i,[x,y] in enumerate(path)]
		self.flux = [y for i,[x,y] in enumerate(path)]
		
		# Move light curve up to avoid negative values
		for i in range(len(self.flux)):
			self.flux[i] += abs(min(self.flux))
	
	def copy(self):
		copy = LightCurveWiener(self.length)
		copy.path = self.path
		for i in range(len(copy.time)):
			copy.time[i] = self.time[i]
			copy.flux[i] = self.flux[i]
		return copy
	
	def regen(self, length):
		self.path = self.Wiener.generate_sample_path(range(length))
		self.time = [x for i,[x,y] in enumerate(self.path)]
		self.flux = [y for i,[x,y] in enumerate(self.path)]
	
	def smooth(self, sigma=1.5):
		import scipy.signal as signal
		
		def gauss_kern(sigma=1.5):
			import numpy as np
			""" Returns a normalized 1D gauss kernel array for convolutions """
			length = int(sigma*10)
			x = np.array(range(length))
			x = x-length/2.
			#sigma = 1.5
			g = np.exp( -( x**2. ) / (2.*sigma**2.) );
			return g / g.sum()
		
		#g = gauss_kern(sigma=len(self.time)/100.)
		g = gauss_kern(sigma)
		
		self.flux = signal.convolve(self.flux,g,mode='same')
	
	def lag_const(self, lag_const):
		for i in range(len(self.time)):
			self.time[i] += lag_const
	
	def lag_luminosity(self, nBins, alpha, beta):
		for i in range(len(self.time)):
			lag = alpha * (self.flux[i])**beta
			self.time[i] += lag
		
		# Since the lags are non linear
		# - luminosities from different epochs will overlap
		# We therefore have to make bins and create a histogram
		
		self.bin(nBins)
	
	def bin(self, nBins):
		binSize = len(self.time)/nBins
		
		binPositions = [0]*nBins
		binValues = [0]*nBins
		
		# Calculate bin positions
		for i in range(nBins):
			binPositions[i] = min(self.time) + i*binSize + 0.5*binSize
		
		# Sum y values for each bin
		for i in range(nBins):
			for j in range(len(self.time)):
				x = self.time[j]
				y = self.flux[j]
				if x > binPositions[i]-0.5*binSize and x < binPositions[i]+0.5*binSize:
					binValues[i] += y
		
		self.time = binPositions
		self.flux = binValues
	
	def scale(self, scale):
		for i in range(len(self.flux)):
			self.flux[i] = self.flux[i] * scale
	
	def plot(self):
		import matplotlib.pyplot as plt
		plt.figure()
		plt.plot(self.time,self.flux)
		#plt.show()
	
	def saveToTxt(self, outFileName='lightcurve.txt'):
		from numpy import savetxt, transpose
		errors = [5.0]*len(self.time)
		savetxt(outFileName, transpose((self.time, self.flux, errors)))


