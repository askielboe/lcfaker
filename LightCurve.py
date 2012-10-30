class LightCurve():
	
	def __init__(self, length=1000):
		
		self.length = length
		
		# # Generate Wiener Path
		# from numpy import asarray
		# import pyprocess as SP
		# # Generate random walk light curve using pyprocesses.py
		# wienerParams = {"mu":0, "sigma":0.1}
		# wienerInitial = {"startTime":0, "startPosition":0, "endTime":self.length, "endPosition": 0 }
		# Wiener = SP.Wiener_process(wienerParams, wienerInitial)
		# self.path = asarray(Wiener.generate_sample_path(range(self.length)))
		
		self.path = self.generateWiener()
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
	
	# def regen(self, length):
	# 	self.path = self.Wiener.generate_sample_path(range(length))
	# 	self.time = [x for i,[x,y] in enumerate(self.path)]
	# 	self.flux = [y for i,[x,y] in enumerate(self.path)]
	
	def smooth(self, sigma=1.5):
		import scipy.signal as signal
		
		length = int(sigma*10)
		
		def gauss_kern(sigma=1.5):
			import numpy as np
			""" Returns a normalized 1D gauss kernel array for convolutions """
			x = np.array(range(length))
			x = x-length/2.
			g = np.exp( -( x**2. ) / (2.*sigma**2.) );
			return g / g.sum()
		
		g = gauss_kern(sigma)
		
		from numpy import asarray
		extraFluxLower = asarray([self.flux[0] for i in range(length)])
		extraFluxUpper = asarray([self.flux[-1] for i in range(length)])
		
		from numpy import concatenate
		extendedFlux = concatenate((extraFluxLower,self.flux,extraFluxUpper))
		
		self.flux = asarray(signal.convolve(extendedFlux,g,mode='same'))
		self.flux = self.flux[length:]
		self.flux = self.flux[:-length]
		
	
	def lag_const(self, lag_const):
		self.time = self.time + lag_const
		print "Running lag_const.."
	
	def lag_luminosity(self, lightCurveCont, c, alpha, beta):
		print "Running lag_luminosity.."
		self.time = self.time + float(c) + float(alpha)*lightCurveCont.flux**float(beta)
		
		# Sort time array
		self.time.sort()
	
	def addNoiseGaussian(self, signalToNoise):
		from random import gauss
		for i in range(len(self.flux)):
			sigma = 1./float(signalToNoise)*self.flux[i]
			self.flux[i] = gauss(self.flux[i],sigma)
	
	def rebin(self, nBins):
		from rebin import rebin
		from numpy import asarray, append
		
		# Since rebin needs bin-edges we have to add an extra point to the time-list
		# Whatever value is in a bin is then identified to lie between the two bin edges
		
		time = self.time
		flux = self.flux
		
		binSize = (max(time)-min(time))/float(nBins)
		tEdges = append(time, max(time)+1.0)
		
		tEdgesNew = asarray([min(time) + i*binSize for i in range(nBins+1)])
		
		self.flux = rebin(tEdges,flux,tEdgesNew)
		
		# To get the correct (x,y) values we have to calulate mid-bin positions based on bin edges:
		from numpy import delete
		tNewMidBin = delete(tEdgesNew,-1)
		tNewMidBin = tNewMidBin + 0.5*binSize
		
		self.time = tNewMidBin
	
	# def bin(self, nBins):
	# 	binSize = (max(self.time)-min(self.time))/float(nBins)
	# 	binPositions = [0]*nBins
	# 	binValues = [0]*nBins
	# 	
	# 	# Calculate bin positions
	# 	for i in range(nBins):
	# 		binPositions[i] = int(min(self.time)) + i*binSize + 0.5*binSize
	# 	
	# 	# Sum y values for each bin
	# 	for i in range(nBins):
	# 		for j in range(len(self.time)):
	# 			x = self.time[j]
	# 			y = self.flux[j]
	# 			if i == nBins-1:
	# 				if x >= binPositions[i]-0.5*binSize and x <= binPositions[i]+0.5*binSize:
	# 					binValues[i] += y
	# 			else:
	# 				if x >= binPositions[i]-0.5*binSize and x < binPositions[i]+0.5*binSize:
	# 					binValues[i] += y
	# 	
	# 	self.time = list(binPositions)
	# 	self.flux = list(binValues)
	
	def scale(self, scale):
		self.flux = self.flux * float(scale)
	
	def plot(self):
		import matplotlib.pyplot as plt
		plt.figure()
		plt.plot(self.time,self.flux)
	
	def getAverageFlux(self):
		return sum(self.flux)/len(self.flux)
	
	def saveToTxt(self, outFileName='lightcurve.txt', signalToNoise=10.):
		from numpy import savetxt, transpose
		errorBars = 1./float(signalToNoise)*self.flux
		savetxt(outFileName, transpose((self.time, self.flux, errorBars)))
