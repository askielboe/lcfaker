def OUprocess(sf, tau, mu, N, tmax):
	import numpy as np
	# Adapted from eq. 5 from MacLeod et al. 2010
	
	dt = tmax / N
	t = np.arange(0.0,tmax,dt)
	
	exp = np.exp(-dt/tau)
	norm = mu * (1-exp)
	var = 0.5 * sf**2.0 * (1-np.exp(-2*dt/tau))
	
	X = np.zeros(N)
	X[0] = np.random.normal(mu,var)
	
	for i in range(N-1):
		E = exp * X[i] + norm
		X[i+1] = np.random.normal(E,var)
	
	return t, X
