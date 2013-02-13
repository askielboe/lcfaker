class LightCurve():
	
	def __init__(self, length=1000):
		
		self.length = length
	
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
	
	def lag_luminosity(self, lightCurveCont, c, alpha, beta):
		print "Running lag_luminosity.."
		
		import numpy as np
		import lib.units as units
		
		# Convert apparant flux to absolute flux before calculating the lag
		absoluteFlux = lightCurveCont.flux - self.mu + self.mag
		print 'mean(absoluteFlux)', np.mean(absoluteFlux)
		
		# Calculate timelag based on radius-luminosity relationship
		timelag = units.r_from_l(units.mag_to_lum5100(absoluteFlux))
		self.time = self.time + timelag
		
		# # #
		# Write some output
		print "Average luminosity: ", np.mean(units.mag_to_lum5100(absoluteFlux))
		
		minLag = units.r_from_l(units.mag_to_lum5100(max(absoluteFlux)))
		maxLag = units.r_from_l(units.mag_to_lum5100(min(absoluteFlux)))
		avgLag = units.r_from_l(units.mag_to_lum5100(np.mean(absoluteFlux)))
		
		print "Minimum lag = ", minLag
		print "Maximum lag = ", maxLag
		print "Lag difference = ", maxLag - minLag
		print "Lag corresponding to average luminosity (Continuum) = ", avgLag
		# # #
		
		# Sort time array
		self.time.sort()
	
	def addNoiseGaussian(self, signalToNoise):
		from random import gauss
		for i in range(len(self.flux)):
			sigma = 1./float(signalToNoise)*self.flux[i]
			self.flux[i] = gauss(self.flux[i],sigma)
	
	def rebin(self, shapeNew):
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
		self.flux = self.flux * float(scale)
	
	def plot(self, units='Jy'):
		import matplotlib.pyplot as plt
		plt.figure()
		plt.title('LightCurve')
		plt.xlabel('time [days]')
		
		if (units == 'Jy'):
			import lib.units as units
			y = units.magi_to_fluxJy(self.flux)
			plt.ylabel('flux [Jy]')
		elif (units == 'Mag'):
			y = self.flux
			plt.ylabel('flux [m_i]')
		else:
			print "ERROR: Unknown unit"
			return
		
		plt.plot(self.time,y)
	
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
	
	def generateOU(self):
		import lcfaker.vendor.pyprocess as SP
		from numpy import asarray
		# Generate random walk light curve using pyprocesses.py
		OUParams = {"theta":1., "mu":5., "sigma":1.25}
		OUInitial = {"startTime":0, "startPosition":0, "endTime":self.length, "endPosition": 0 }
		OU = SP.OU_process(OUParams, OUInitial)
		self.path = asarray(OU.generate_sample_path(range(self.length)))

class LightCurveMacLeod(LightCurve):
	import numpy as np
	
	time = np.array([])
	flux = np.array([])
	
	def __init__(self, \
		nDays = 1000, maxDays = 1000, \
		magApparant = 15, magAbsolute = -24.0, \
		mass = 1e9, lambdarf = 5100.0, z = 0.0):
		
		self.nDays = nDays
		self.maxDays = maxDays
		self.magApparant = magApparant
		self.magAbsolute = magAbsolute
		self.mass = mass
		self.lambdarf = lambdarf
		self.z = z
		
		# Get random walk parameters from physical parameters
		self.sf, self.tau = self.calcSFandTau()
		
		# Generate lightcurve based on damped random walk using equations from MacLeod et al.
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
			+ C * (self.magAbsolute + 23.) \
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
			+ C * (self.magAbsolute + 23.) \
			+ D * np.log(self.mass/1e9) \
			+ E * np.log(1 + self.z)
	
		tau = np.exp(logtau)
	
		return (sf,tau)
	
	def generateMacLeod(self):
		from lib.ouprocess import OUprocess
		self.time, self.flux = OUprocess(self.nDays, self.maxDays, self.magApparant, self.sf, self.tau)

