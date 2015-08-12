def magi_to_lum5100(M):
	# Calculate luminosity L(5100 A)
	Msun = 4.83 # Sun absolute magnitude
	Lsun = 3.846e26
	#M = -2.5*np.log(L/Lsun) + Msun
	#L = 10.0**(0.4*(Msun-M))*Lsun
	L = 10.0**(0.4*(Msun-M))*Lsun*5100.0
	return L

def r_from_l(L):
	import numpy as np
	# Radius-luminosity relation from Bentz et al.
	K = -21.3 # + 2.9 - 2.8
	alpha = 0.519 # + 0.063 - 0.066
	R = 10.0**(K) * 10.0**(alpha*np.log10(L))
	return R

def magi_to_flux(mag):
	# Based on: https://en.wikipedia.org/wiki/Apparent_magnitude
	F_i_0 = 4.76
	F = F_i_0 * 10.0**(-0.4*mag)
	return F*1e-20

def magr_to_flux(mag):
	# Based on: https://en.wikipedia.org/wiki/Apparent_magnitude
	F_r_0 = 4.49
	F = F_r_0 * 10.0**(-0.4*mag)
	return F*1e-20

def magi_to_fluxJy(mag):
	# Based on: https://en.wikipedia.org/wiki/Apparent_magnitude
	F_i_0 = 4760.0
	F = F_i_0 * 10.0**(-0.4*mag)
	return F

# def lag_from_mag(mag):
# 	# Based on

def sdss_r_to_i_absolute(mag_r):
	# Straight line fit to SDSS data
	# Converted to absolute mags using LambdaCDM
	a = 0.95880946699563219
	b = -1.1830349097544115
	return a*mag_r + b

	299792458.0/7900e-10 * 4.76