def mag_to_lum5100(M):
	# Calculate luminosity L(5100 A)
	Msun = 4.83 # Sun absolute magnitude
	Lsun = 3.846e26
	#M = -2.5*np.log(L/Lsun) + Msun
	#L = 10.0**(0.4*(Msun-M))*Lsun
	FudgeFactor = 1000000
	L = 10.0**(0.4*(Msun-M))*Lsun*FudgeFactor
	return L

def r_from_l(L):
	import numpy as np
	# Radius-luminosity relation from Bentz et al.
	K = -21.3 # + 2.9 - 2.8
	alpha = 0.519 # + 0.063 - 0.066
	R = 10.0**(K) * 10.0**(alpha*np.log10(L))
	return R
