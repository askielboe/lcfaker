class LightCurve():
	
	def __init__(self, length=1000):
		
		self.length = length
	
	# def regen(self, length):
	# 	self.path = self.Wiener.generate_sample_path(range(length))
	# 	self.time = [x for i,[x,y] in enumerate(self.path)]
	# 	self.flux = [y for i,[x,y] in enumerate(self.path)]
	
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
		self.time = self.time + lag_const
		print "Running lag_const.."
	
	def lag_luminosity(self, lightCurveCont, c, alpha, beta):
		import lib.physics as phys
		print "Running lag_luminosity.."
		
		self.time = self.time + phys.radius_from_luminosity_relation(phys.mag_to_lum5100(lightCurveCont.flux))
		
		#self.time = self.time + float(c) + float(alpha)*lightCurveCont.flux**float(beta)
		
		# Sort time array
		self.time.sort()
	
	def addNoiseGaussian(self, signalToNoise):
		from random import gauss
		for i in range(len(self.flux)):
			sigma = 1./float(signalToNoise)*self.flux[i]
			self.flux[i] = gauss(self.flux[i],sigma)
	
	def rebin(self, nBins):
		from lcfaker.vendor.rebin import rebin
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

class LightCurveSP(LightCurve):
	def __init__(self, length=1000):
		self.length = length
		self.path = self.generateWiener()
		self.time = self.path[:,0]
		self.flux = self.path[:,1]
	
	def generateWiener(self):
		import lcfaker.vendor.pyprocess as SP
		from numpy import asarray
		# Generate random walk light curve using pyprocesses.py
		wienerParams = {"mu":0, "sigma":0.1}
		wienerInitial = {"startTime":0, "startPosition":0, "endTime":self.length, "endPosition": 0 }
		Wiener = SP.Wiener_process(wienerParams, wienerInitial)
		return asarray(Wiener.generate_sample_path(range(self.length)))
		#self.time = [x for i,[x,y] in enumerate(path)]
		#self.flux = [y for i,[x,y] in enumerate(path)]
	
	def generateOU(self):
		import lcfaker.vendor.pyprocess as SP
		from numpy import asarray
		# Generate random walk light curve using pyprocesses.py
		OUParams = {"theta":1., "mu":5., "sigma":1.25}
		OUInitial = {"startTime":0, "startPosition":0, "endTime":self.length, "endPosition": 0 }
		OU = SP.OU_process(OUParams, OUInitial)
		self.path = asarray(OU.generate_sample_path(range(self.length)))
		#self.time = [x for i,[x,y] in enumerate(path)]
		#self.flux = [y for i,[x,y] in enumerate(path)]

class LightCurveMacLeod(LightCurve):
	path = [[],[]]
	time = []
	flux = []
	
	def __init__(self, mu = -23.0, mag = -23.0, mass = 1e9, lambdarf = 5100.0, z = 0.0):
		#self.length = length
		#self.path = self.generateMacLeod(-23.0)
		#self.time = []
		#self.flux = []
		
		self.N = 1000
		self.tmax = 3650.
		
		self.mu = mu
		self.mag = mag
		self.mass = mass
		self.lambdarf = lambdarf
		self.z = z
		
		self.sf, self.tau = self.calcSFandTau()
		
		self.generateMacLeod()
	
	def calcSFandTau(self):
		import numpy as np
		# Relation from MacLeod et al. 2010 eq. 7
		# Parameters from MacLeod et al. 2010 Table 1
		A = -0.51
		B = -0.479
		C = 0.131
		D = 0.18
		E = 0.0
	
		logsf = A \
			+ B * np.log(self.lambdarf/4000.) \
			+ C * (self.mag + 23.) \
			+ D * np.log(self.mass/1e9) \
			+ E * np.log(1 + self.z)
	
		sf = np.exp(logsf)
	
		A = 2.4
		B = 0.17
		C = 0.03
		D = 0.21
		E = 0.0
	
		logtau = A \
			+ B * np.log(self.lambdarf/4000.) \
			+ C * (self.mag + 23) \
			+ D * np.log(self.mass/1e9) \
			+ E * np.log(1 + self.z)
	
		tau = np.exp(logtau)
	
		return sf,tau
	
	def generateMacLeod(self):
		from lib.ouprocess import OUprocess
		self.time, self.flux = OUprocess(self.sf, self.tau, self.mu, self.N, self.tmax)
		self.path = [self.time,self.flux]
