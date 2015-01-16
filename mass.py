""" File to calculate total mass of galaxy from (u-r) colours and black hole masses from FWHM of Halpha broad lines in galaxy spectra."""

import numpy as N
import pylab as P
import pyfits as F

from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=71.0, Om0 = 0.26)
H0 = cosmo.H0.value

msol = 1.989E30
mpc = 3.08567758E16 * 1E6 # mpc in m
c = 299792.458 #speed of light in km/s



def calc_ur(psf_u, psfErr_u, psf_r, psfErr_r, petro_u, petroErr_u, petro_r, petroErr_r):
	u = petro_u - psf_u
	Err_u = (psfErr_u **2 + petroErr_u**2) ** 0.5
	r = petro_r - psf_r
	Err_r = (psfErr_r **2 + petroErr_r**2) ** 0.5
	ur = u - r
	Err_ur = (Err_u **2 + Err_r**2) ** 0.5
	return ur, Err_ur, petro_r

def bhmass(z, flux, fwhm):
    #d = cosmo.comoving_distance(z).value * mpc * 100 # distance in centimetres
    d = (z*c)/H0 * mpc * 100 # distance in centimeters
    print 'distance', d
    #ew_v = (ew/6562.8) * c #ew_v in km/s
    fwhm_v = (fwhm/6562.8) * c #fwhm in km/s
    print 'fwhm', fwhm
    lum = flux*1E-17*6562.8*(4*N.pi*(d**2)) # luminosity in erg/s but flux in erg/cm2/s/A - so integrate over all wavelengths - i.e. multiply by the wavelength
    print 'luminosity', lum
    mbh = 3E6 * ((lum/1E42)**0.45) * ((fwhm_v/1E3)**2.06)
    return N.log10(mbh)

 
def calc_fwhm(wave, emission):
    max = N.max(emission)
    i = N.argmax(emission)
    hm = max/2
    idx = N.abs(emission-hm)
    fidx = N.argmin(idx[:i])
    sidx = N.argmin(idx[i:]) + i
    fwhm = wave[sidx] - wave[fidx]
    return fidx, sidx, fwhm

def gauss(a, u, s, wave):
        return a * N.exp(- (((wave - u)**2) /(s**2)))

sdss = N.genfromtxt('/Volumes/Data/smethurst/observing/la_palma/int_reduc/result.csv', delimiter=',', skip_header=True)

sourceList = ['Q078017', 'Q507953', 'Q464967', 'Q292635', 'Q186225']
sourceID = ['1237679457600078017', '1237666245208507953', '1237668631602464967', '1237662981045292635', '1237671688011186225']
target_flux = [11.34, 50.72, 36.55, 70.92, 672.08]
target_fwhm = [0.3, 2.92, 9.51, 6.61, 24.90]
emission_fwhm = [0.3, 2.92, 9.51, 6.61, 24.90]
target_I_v0 = N.zeros(len(sourceList))

sdss = N.append(sdss, N.zeros((5,6)), axis=1)

z = N.zeros_like(sourceList)
bh_masses = N.zeros(len(sourceList))

for n in range(len(z)):
    s = F.open('/Volumes/Data/smethurst/observing/la_palma/data_may14/'+sourceList[n]+'_extract_agn.fit')
    d = s[0].data
    h = s[0].header
    #l = s[0].header['crval1'] + s[0].header['cd1_1']*(N.arange(s[0].header['naxis1']))
    w = N.linspace(0-h['CRPIX2'], d.shape[0]-h['CRPIX2'], d.shape[0])
    l = w*h['cd2_2'] + h['crval2']
    if sourceList[n] == 'Q464967':
        ml = l[N.argmax(d[100:-100])+100]
        target_I_v0[n] = N.max(d[100:-100]) / 1E-17
    if sourceList[n] == 'Q078017':
        ml = l[N.argmax(d[:100])]
        target_I_v0[n] = N.max(d[:100]) / 1E-17
    else:
        ml = l[N.argmax(d)]
        target_I_v0[n] = N.max(d) / 1E-17
    z[n] = (ml/6562.8) - 1
    print 'lambda value', ml
    print 'z = ',z[n]
    sdss[n,17] = z[n]   

for n in range(len(sourceList)):
    ur, Err_ur, r = sdss[n, -5:-2] = calc_ur(sdss[n,5], sdss[n,6], sdss[n, 7], sdss[n,8], sdss[n,13], sdss[n,14], sdss[n,15], sdss[n,16])
    print 'ur ', ur
    if ur <= 2.1:
    	log_m_l = -0.95 + 0.56 * ur
    else:
        log_m_l = -0.16 + 0.18 * ur
    ld = cosmo.luminosity_distance(sdss[n,17]).value * 1E6
    Mr = r - 5 * (N.log10(ld) - 1)
    print 'Mr ', Mr
    m_msun = sdss[n, -2] = ((4.62 - Mr)/2.5) + log_m_l
    bf = F.open('/Volumes/Data/smethurst/idl/gandalf_release_v1.5/SDSS_example/'+sourceList[n]+'_extract_agn_deredshift_rebin_header_units_GANDALF_fits.fits')
    hdr = bf[0].header
    emission = bf[2].data
    lam = hdr['crval1'] + hdr['cd1_1']*(N.arange(hdr['naxis1'] - hdr['crpix1']))
    wave = 10**lam
    if sourceList[n] == 'Q078017':
        target_fwhm[n] = calc_fwhm(wave, emission)
    if sourceList[n] = 'Q507953'
        target_fwhm[n] = calc_fwhm(wave, emission)
    else:
        bf = N.load('best_fit_'+source[n]+'.npy')
        broad = gauss(bf[3][0], bf[4][0], bf[5][0], wave)
        target_fwhm[n] = calc_fwhm(wave, broad)
    print target_fwhm
    emission_fwhm[n] = calc_fwhm(wave, emission)
    print emission_fwhm
    bh_masses[n] = sdss[n, -1] = bhmass(sdss[n,17], target_I_v0[n], target_fwhm[n])

print 'u-r', sdss[:,-5]

P.figure()
P.scatter(sdss[:,-2], sdss[:,-1])
P.xlabel(r'$log_{10}(M_{*}/M_{\odot})$')
P.ylabel(r'$log_{10}(M_{BH}/M_{\odot})$')
P.show()



