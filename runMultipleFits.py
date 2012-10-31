### Load dependencies
import numpy as np

import FakeLight
from calcPosteriors import calcPosteriors

from javelin.zylc import get_data
from javelin.lcmodel import Cont_Model, Rmap_Model


### Define running parameters

nRuns = 1

signalToNoise = 10.

observeIntervalCenters = [250,400,500,650,800]
observeIntervalWidths = [45,55,30,60,65]

cList = [50.]*nRuns
alphaList = [float(i+1.) for i in range(nRuns)]
betaList = [0.519]*nRuns

### Define result parameters

meanList = [0.]*nRuns
stdList = [0.]*nRuns

### Functions

def runFakeLight(c, alpha, beta):
	# Generate light curve
	lc = FakeLight.FakeLight(1000)
	
	lc.reprocess(200, 1.0, 1.5, float(c), float(alpha), float(beta))
	
	lc.lightCurveCont.addNoiseGaussian(signalToNoise)
	lc.lightCurveLine.addNoiseGaussian(signalToNoise)
	
	lc.observeIntervals(observeIntervalCenters,observeIntervalWidths)
	
	lc.saveToTxt()

def runJavelin():
	# Run JAVELIN
	c = get_data(["lc_cont.txt"])
	cmod = Cont_Model(c)
	cmod.do_mcmc()
	
	cy = get_data(["lc_cont.txt","lc_line.txt"], names=["Continuum", "Line"])
	cymod = Rmap_Model(cy)
	cymod.do_mcmc(conthpd=cmod.hpd)
	
	# Save chains
	np.savetxt('cymod_flatchain.txt', cymod.flatchain)

### Run FakeLight and JAVELIN
for i in range(nRuns):
	### Run with constant lag
	runFakeLight(cList[i], alphaList[i], betaList[i])
	runJavelin()
	
	# Get and append t-lag posterior
	(mean,std) = calcPosteriors('cymod_flatchain.txt', 0., 100.)
	meanList[i] = mean
	stdList[i] = std
	
	# ### Run with L dependent lag
	# runFakeLight(20., 11.0, 0.519)
	# runJavelin()
	# 
	# # Get and append t-lag posterior
	# (mean,std) = calcPosteriors('cymod_flatchain.txt', 0., 100.)
	# meanList[i] = mean
	# stdList[i] = std

# print(cList)
# print(alphaList)
# print(betaList)
# print(meanList)
# print(stdList)

outArray = np.transpose(np.array([cList,alphaList,betaList,meanList,stdList]))
np.savetxt('runMultipleFitsOutput.txt', outArray)
