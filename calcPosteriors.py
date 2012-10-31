def calcPosteriors(filename, xMin, xMax):

	import numpy as np
	
	data = np.loadtxt(filename)
	
	tLag = data[:,2]
	
	# Remove everything outside selected range
	tLag.sort()
	
	xMin = 0.
	xMax = 100.
	
	iMin = min(np.where(tLag >= xMin)[0])
	iMax = max(np.where(tLag <= xMax)[0])
	
	tLagMasked = tLag[iMin:iMax]
	
	mean = tLagMasked.mean()
	std = tLagMasked.std()
	
	return (mean,std)
	
	# import scipy.stats as stats
	# gaussian = stats.norm
	# 
	# mean, std = gaussian.fit(tLagMasked)
	# print(mean, std)
