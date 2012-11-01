### Load dependencies
import numpy as np

import FakeLight
from calcPosteriors import calcPosteriors

from javelin.zylc import get_data
from javelin.lcmodel import Cont_Model, Rmap_Model

### Define running parameters
nRuns = 11 # Range of different parameters
nSamples = 10 # Number of runs with same parameters to deterine scatter

signalToNoiseList = [3.]*nRuns
#signalToNoiseList = [float(i+1.) for i in range(nRuns)]

cList = [50.]*nRuns

#alphaList = [0.]*nRuns
alphaList = [float(i) for i in range(nRuns)]

betaList = [0.519]*nRuns

observeIntervalCenters = [250,400,500,650,800]
observeIntervalWidths = [45,55,30,60,65]

outFileName = "alpha_"+str(min(alpha))+"_"+str(max(alpha))+"_SN_"+str(min(signalToNoiseList))+"_"+str(max(signalToNoiseList))+".txt"

### Functions

def runFakeLight(i, c, alpha, beta, signalToNoise):
	# Generate light curve
	lc = FakeLight.FakeLight(1000)
	
	(minLag,maxLag,avgLag) = lc.reprocess(200, 1.0, 1.5, float(c), float(alpha), float(beta))
	
	lc.lightCurveCont.addNoiseGaussian(signalToNoise)
	lc.lightCurveLine.addNoiseGaussian(signalToNoise)
	
	lc.observeIntervals(observeIntervalCenters,observeIntervalWidths)
	
	lc.saveToTxt("lc_cont_"+str(i)+".txt","lc_line_"+str(i)+".txt")
	
	return (minLag, maxLag, avgLag)

def runJavelin(i):
	# Run JAVELIN
	c = get_data(["lc_cont_"+str(i)+".txt"])
	cmod = Cont_Model(c)
	cmod.do_mcmc()
	
	cy = get_data(["lc_cont_"+str(i)+".txt","lc_line_"+str(i)+".txt"], names=["Continuum", "Line"])
	cymod = Rmap_Model(cy)
	cymod.do_mcmc(conthpd=cmod.hpd)
	
	# Save chains
	np.savetxt("cymod_flatchain_"+str(i)+".txt", cymod.flatchain)

def runFakeAndJavelin(i):
	
	meanSamples = [0.]*nSamples
	stdSamples = [0.]*nSamples
	
	minLagSamples = [0.]*nSamples
	maxLagSamples = [0.]*nSamples
	avgLagSamples = [0.]*nSamples
	
	for j in range(nSamples):
		(minLag, maxLag, avgLag) = runFakeLight(i, cList[i], alphaList[i], betaList[i], signalToNoiseList[i])
		
		minLagSamples[j] = minLag
		maxLagSamples[j] = maxLag
		avgLagSamples[j] = avgLag
		
		runJavelin(i)
		
		# Get and append t-lag posterior
		(mean,std) = calcPosteriors("cymod_flatchain_"+str(i)+".txt", 0., 100.)
		
		meanSamples[j] = mean
		stdSamples[j] = std
	
	meanMean = np.mean(meanSamples)
	meanStd = np.std(meanSamples)
	
	stdMean = np.mean(stdSamples)
	stdStd = np.std(stdSamples)
	
	minLagMean = np.mean(minLagSamples)
	maxLagMean = np.mean(maxLagSamples)
	avgLagMean = np.mean(avgLagSamples)
	
	minLagStd = np.std(minLagSamples)
	maxLagStd = np.std(maxLagSamples)
	avgLagStd = np.std(avgLagSamples)
	
	return (
		meanMean,
		meanStd,
		stdMean,
		stdStd,
		minLagMean,
		maxLagMean,
		avgLagMean,
		minLagStd,
		maxLagStd,
		avgLagStd
		)

### Run FakeLight and JAVELIN

### Define result parameters

meanMeanList = [0.]*nRuns
meanStdList = [0.]*nRuns

stdMeanList = [0.]*nRuns
stdStdList = [0.]*nRuns

minLagMeanList = [0.]*nRuns
maxLagMeanList = [0.]*nRuns
avgLagMeanList = [0.]*nRuns

minLagStdList = [0.]*nRuns
maxLagStdList = [0.]*nRuns
avgLagStdList = [0.]*nRuns

## Using multiprocessing
from multiprocessing import Pool
pool = Pool(nRuns)

(
	meanMeanList,
	meanStdList,
	stdMeanList,
	stdStdList,
	minLagMeanList,
	maxLagMeanList,
	avgLagMeanList,
	minLagStdList,
	maxLagStdList,
	avgLagStdList
	) = zip(*pool.map(runFakeAndJavelin, range(nRuns)))



# for i in range(nRuns):
# 	for j in range(nSamples):
# 		### Run with constant lag
# 		(minLag, maxLag, avgLag) = runFakeLight(cList[i], alphaList[i], betaList[i], signalToNoiseList[i])
# 		
# 		minLagSamples[j] = minLag
# 		maxLagSamples[j] = maxLag
# 		avgLagSamples[j] = avgLag
# 		
# 		runJavelin()
# 		
# 		# Get and append t-lag posterior
# 		(mean,std) = calcPosteriors("cymod_flatchain_"+str(i)+".txt", 0., 100.)
# 		meanSamples[j] = mean
# 		stdSamples[j] = std
# 	
# 	meanMeanList[i] = np.mean(meanSamples)
# 	meanStdList[i] = np.std(meanSamples)
# 	
# 	stdMeanList[i] = np.mean(stdSamples)
# 	stdStdList[i] = np.std(stdSamples)
# 	
# 	minLagMeanList[i] = np.mean(minLagSamples)
# 	maxLagMeanList[i] = np.mean(maxLagSamples)
# 	avgLagMeanList[i] = np.mean(avgLagSamples)
# 	
# 	minLagStdList[i] = np.std(minLagSamples)
# 	maxLagStdList[i] = np.std(maxLagSamples)
# 	avgLagStdList[i] = np.std(avgLagSamples)
# 	
# 	meanSamples = [0.]*nSamples
# 	stdSamples = [0.]*nSamples
# 	
# 	minLagSamples = [0.]*nSamples
# 	maxLagSamples = [0.]*nSamples
# 	avgLagSamples = [0.]*nSamples

outArray = np.transpose(np.array([
	cList,
	alphaList,
	betaList,
	minLagMeanList,
	minLagStdList,
	maxLagMeanList,
	maxLagStdList,
	avgLagMeanList,
	avgLagStdList,
	signalToNoiseList,
	meanMeanList,
	meanStdList,
	stdMeanList,
	stdStdList]))

np.savetxt(outFileName, outArray)
