""" File to calculate total mass of galaxy from (u-r) colours and black hole masses from FWHM of Halpha broad lines in galaxy spectra."""

import numpy as N
import pylab as P
import pyfits as F

# class Galaxy():
# 	pass

# Q078017 = Galaxy()
# Q078017.id = '1237679457600078017'
# Q078017.ra = 337.581875
# Q078017.dec = 16.520222 

# Q508953 = Galaxy()
# Q508953.id = '1237666245208507953'
# Q508953.ra = 308.663542
# Q508953.dec = 14.569917

# Q464967 = Galaxy()
# Q464967.id = '1237668631602464967'
# Q464967.ra = 262.202000
# Q464967.dec = 0.798333

# Q292365 = Galaxy()
# Q292365.id = '1237662981045292635'
# Q292365.ra = 290.711458
# Q292365.dec = -5.883028

# Q186225 = Galaxy()
# Q186225.id = '1237671688011186225'
# Q186225.ra = 178.409375
# Q186225.dec = 72.864944 

sdss = N.genfromtxt('/Volumes/Data/smethurst/observing/la_palma/int_reduc/result.csv', delimiter=',', names=['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm'])

sourceList = [Q078017, Q508953, Q464967, Q292365, Q186225]
sourceID = ['1237679457600078017', '1237666245208507953', '1237668631602464967', '1237662981045292635', '1237671688011186225'] 

def calc_ur(psf_u, psfErr_u, psf_r, psfErr_r, petro_u, petroErr_u, petro_r, petroErr_r):
	u = petro_u - psf_u
	Err_u = (psfErr_u **2 + petroErr_u**2) ** 0.5
	r = petro_r - psf_r
	Err_r = (psfErr_r **2 + petroErr_r**2) ** 0.5
	ur = u - r
	Err_ur = (Err_u **2 + Err_r**2) ** 0.5
	return ur, Err_ur

sdss = N.append(sdss, N.zeros((5,2)), axis=1)
for n in range(len(sourceList)):
	sdss[n, -2:] = calc_ur(sdss[n,5], sdss[n,6], sdss[n, 7], sdss[n,8], sdss[n,9], sdss[n,10], sdss[n,11], sdss[n,12])

