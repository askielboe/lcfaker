def OUprocess(N, Nmax, mu, sf, tau):
	import numpy as np
	# Adapted from eq. 5 from MacLeod et al. 2010

	dt = Nmax / N
	t = np.arange(0.0,Nmax,dt)

	exp = np.exp(-dt/tau)
	norm = mu * (1.0-exp)
	var = 0.5 * sf**2.0 * (1.0-np.exp(-2.0*dt/tau))
	std = np.sqrt(var)

	X = np.zeros(N)
	X[0] = np.random.normal(mu,std)

	for i in range(N-1):
		E = exp * X[i] + norm
		X[i+1] = np.random.normal(E,std)

	return (t, X)
