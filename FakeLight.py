class LightContAndLine():
	
	timeCont = []
	fluxCont = []
	timeLine = []
	fluxLine = []
	
	def __init__(self, length):
		
		# Input parameters
		self.length = length
		# self.nBins = nBins
		# self.scale = scale
		# self.sigma = sigma
		# self.alpha = alpha
		# self.beta = beta
		
		# Generate a lightcurve
		self.lightCurve = LightCurve(self.length)
		
		# Now done when defining the light curve
		#lightCurve.generateWiener()
		#lightCurve.generateOU()
		
		# Normalize
		
		# for i in range(len(self.lightCurve.flux)):
			#self.lightCurve.flux[i] = self.lightCurve.flux[i] / max(self.lightCurve.flux) * 5. + 5.
		
		# Normalize using numpy array
		self.lightCurve.path[:,1] = self.lightCurve.path[:,1] / max(abs(self.lightCurve.path[:,1])) + 1.
		
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
		
		self.lightCurveCont.time = timeObserved
		self.lightCurveLine.time = timeObserved
		
		self.lightCurveCont.flux = fluxContObserved
		self.lightCurveLine.flux = fluxLineObserved
		
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

class LightCurve():
	#from numpy import array
	
	#time = array(0)
	#flux = array(0)
	
	def __init__(self, length):
		
		self.length = length
		
		# Generate Wiener Path
		from numpy import asarray
		import pyprocess as SP
		# Generate random walk light curve using pyprocesses.py
		wienerParams = {"mu":0, "sigma":0.1}
		wienerInitial = {"startTime":0, "startPosition":0, "endTime":self.length, "endPosition": 0 }
		Wiener = SP.Wiener_process(wienerParams, wienerInitial)
		self.path = asarray(Wiener.generate_sample_path(range(self.length)))
		self.time = self.path[:,0]
		self.flux = self.path[:,1]
	
	def generateWiener(self):
		import pyprocess as SP
		from numpy import asarray
		# Generate random walk light curve using pyprocesses.py
		wienerParams = {"mu":0, "sigma":0.1}
		wienerInitial = {"startTime":0, "startPosition":0, "endTime":self.length, "endPosition": 0 }
		Wiener = SP.Wiener_process(wienerParams, wienerInitial)
		return asarray(Wiener.generate_sample_path(range(self.length)))
		#self.time = [x for i,[x,y] in enumerate(path)]
		#self.flux = [y for i,[x,y] in enumerate(path)]
		
		# Move light curve up to avoid negative values
		for i in range(len(self.flux)):
			self.flux[i] += abs(min(self.flux))
	
	def generateOU(self):
		import pyprocess as SP
		from numpy import asarray
		# Generate random walk light curve using pyprocesses.py
		OUParams = {"theta":1., "mu":5., "sigma":1.25}
		OUInitial = {"startTime":0, "startPosition":0, "endTime":self.length, "endPosition": 0 }
		OU = SP.OU_process(OUParams, OUInitial)
		self.path = asarray(OU.generate_sample_path(range(self.length)))
		#self.time = [x for i,[x,y] in enumerate(path)]
		#self.flux = [y for i,[x,y] in enumerate(path)]
		
		# Move light curve up to avoid negative values
		#for i in range(len(self.flux)):
		#	self.flux[i] += abs(min(self.flux))
	
	def copy(self):
		copy = LightCurve(self.length)
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
		
		length = int(sigma*10)
		
		def gauss_kern(sigma=1.5):
			import numpy as np
			""" Returns a normalized 1D gauss kernel array for convolutions """
			x = np.array(range(length))
			x = x-length/2.
			#sigma = 1.5
			g = np.exp( -( x**2. ) / (2.*sigma**2.) );
			return g / g.sum()
		
		#g = gauss_kern(sigma=len(self.time)/100.)
		g = gauss_kern(sigma)
		
		from numpy import asarray
		extraFluxLower = asarray([self.flux[0] for i in range(length)])
		extraFluxUpper = asarray([self.flux[-1] for i in range(length)])
		
		from numpy import concatenate
		extendedFlux = concatenate((extraFluxLower,self.flux,extraFluxUpper))
		
		self.flux = signal.convolve(extendedFlux,g,mode='same')
		
		#self.flux = signal.convolve(extraFluxLower+self.flux+extraFluxUpper,g,mode='same')
		
		self.flux = self.flux[length:]
		self.flux = self.flux[:-length]
		
		#self.flux = signal.convolve(self.flux,g,mode='same')
	
	def lag_const(self, lag_const):
		#for i in range(len(self.time)):
			#self.time[i] += lag_const
		
		# Using numpy array
		#self.time += lag_const
		self.time = self.time + lag_const
		print "Running lag_const.."
	
	def lag_luminosity(self, lightCurveCont, c, alpha, beta):
		# for i in range(len(self.time)):
		# 	lag = c + alpha*(lightCurveCont.flux[i])**beta
		# 	self.time[i] += lag
		# Using ndarray
		print "Running lag_luminosity.."
		self.time = self.time + float(c) + float(alpha)*lightCurveCont.flux**float(beta)
		
		# Sort time array
		self.time.sort()
	
	def rebin(self, nBins):
		from rebin import rebin
		from numpy import asarray, append
		
		# Convert lists to arrays
		# Since rebin needs bin-edges we have to add an extra point to the time-list
		# Whatever value is in a bin is then identified to lie between the two bin edges
		
		time = self.time
		flux = self.flux
		
		binSize = (max(time)-min(time))/float(nBins)
		#tEdges = concatenate((time,asarray(max(time)+1.0)))
		tEdges = append(time, max(time)+1.0)
		
		tEdgesNew = asarray([min(time) + i*binSize for i in range(nBins+1)])
		
		self.flux = rebin(tEdges,flux,tEdgesNew)
		
		# To get the correct (x,y) values we have to calulate mid-bin positions based on bin edges:
		from numpy import delete
		tNewMidBin = delete(tEdgesNew,-1)
		tNewMidBin = tNewMidBin + 0.5*binSize
		
		self.time = tNewMidBin
	
	def bin(self, nBins):
		binSize = (max(self.time)-min(self.time))/float(nBins)
		binPositions = [0]*nBins
		binValues = [0]*nBins
		
		# Calculate bin positions
		for i in range(nBins):
			binPositions[i] = int(min(self.time)) + i*binSize + 0.5*binSize
		
		# Sum y values for each bin
		for i in range(nBins):
			for j in range(len(self.time)):
				x = self.time[j]
				y = self.flux[j]
				if i == nBins-1:
					if x >= binPositions[i]-0.5*binSize and x <= binPositions[i]+0.5*binSize:
						binValues[i] += y
				else:
					if x >= binPositions[i]-0.5*binSize and x < binPositions[i]+0.5*binSize:
						binValues[i] += y
		
		self.time = list(binPositions)
		self.flux = list(binValues)
	
	def scale(self, scale):
		#for i in range(len(self.flux)):
			#self.flux[i] = self.flux[i] * scale
		# Using ndarray
		self.flux = self.flux * float(scale)
	
	def plot(self):
		import matplotlib.pyplot as plt
		plt.figure()
		plt.plot(self.time,self.flux)
		#plt.show()
	
	def getAverageFlux(self):
		return sum(self.flux)/len(self.flux)
	
	def saveToTxt(self, outFileName='lightcurve.txt'):
		from numpy import savetxt, transpose
		errors = [max(self.flux)/10.]*len(self.time)
		savetxt(outFileName, transpose((self.time, self.flux, errors)))


