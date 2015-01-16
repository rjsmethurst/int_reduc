"""File to use emcee to fit a double gaussian with a narrow and broad peak to emission data from GANDALF."""

import numpy as N
import pylab as P
import pyfits as F
import emcee 
import triangle
import time

def walker_plot(samples, nwalkers, limit):
        s = samples.reshape(nwalkers, -1, 6)
        s = s[:,:limit, :]
        fig = P.figure(figsize=(8,15))
        ax1 = P.subplot(6,1,1)
        ax2 = P.subplot(6,1,2)
        ax3 = P.subplot(6,1,3)
        ax4 = P.subplot(6,1,4)
        ax5 = P.subplot(6,1,5)
        ax6 = P.subplot(6,1,6)
        for n in range(len(s)):
            ax1.plot(s[n,:,0], 'k')
            ax2.plot(s[n,:,1], 'k')
            ax3.plot(s[n,:,2], 'k')
            ax4.plot(s[n,:,3], 'k')
            ax5.plot(s[n,:,4], 'k')
            ax6.plot(s[n,:,5], 'k')
        ax1.tick_params(axis='x', labelbottom='off')
        ax2.tick_params(axis='x', labelbottom='off')
        ax3.tick_params(axis='x', labelbottom='off')
        ax4.tick_params(axis='x', labelbottom='off')
        ax5.tick_params(axis='x', labelbottom='off')
        P.subplots_adjust(hspace=0.1)
        save_fig = 'walkers_steps_'+str(time.strftime('%H_%M_%d_%m_%y'))+'.pdf'
        fig.savefig(save_fig)
        return fig

def model(theta, wave):
	a, u, s, a1, u1, s1 = theta
	return a * N.exp(- (((wave - u)**2) /(s**2))) + a1 * N.exp(- (((wave - u1)**2) /(s1**2)))

def lnlike(theta, wave, flux, fluxerr):
	pred_f = model(theta, wave)
	chi =  -0.5*N.log(2*N.pi*fluxerr**2)-0.5*((flux-pred_f)**2/fluxerr**2)
    	return N.sum(chi)

def lnprior(theta, wave, flux):
	a, u, s, a1, u1, s1 = theta
	if N.min(flux) < a < N.max(flux)+1 and N.min(flux) < a1 < N.max(flux)+1 and wave[0] < u < wave[-1] and wave[0] < u1 < wave[-1] and 0 < s < 1000 and 0 < s1 < 1000 and s < s1:
		return 0.0
	else:
		return -N.inf

def lnprob(theta, wave, flux, fluxerr):
	lp = lnprior(theta, wave, flux)
    	if not N.isfinite(lp):
        	return -N.inf
    	return lp + lnlike(theta, wave, flux, fluxerr)

ndim = 6
nwalkers = 100
nsteps = 2000
burnin = 5000
start = [6, 6570, 2, 6, 6570, 30]

#source = ['Q078017', 'Q507953', 'Q464967', 'Q292635', 'Q186225']
source =['Q078017', 'Q507953']

for n in range(len(source)):
	bf = F.open('/Volumes/Data/smethurst/idl/gandalf_release_v1.5/SDSS_example/'+source[n]+'_extract_agn_deredshift_rebin_header_units_GANDALF_fits.fits')
    	hdr = bf[0].header
    	flux = bf[2].data
    	fluxerr = N.random.randn(len(flux))
    	lam = hdr['crval1'] + hdr['cd1_1']*(N.arange(hdr['naxis1'] - hdr['crpix1']))
    	wave = 10**lam
	p0 = [start + 1e-4*N.random.randn(ndim) for i in range(nwalkers)]
	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(wave, flux, fluxerr))
	# burn in run 
	pos, prob, state = sampler.run_mcmc(p0, burnin)
	samples = sampler.chain[:,:,:].reshape((-1,ndim))
	N.save('samples_burn_in_'+source[n]+'.npy', samples)
	walker_plot(samples, nwalkers, burnin)
	sampler.reset()
	print 'RESET', pos
	# main sampler run 
	sampler.run_mcmc(pos, nsteps)
	samples = sampler.chain[:,:,:].reshape((-1,ndim))
	N.save('samples_'+source[n]+'.npy', samples)
	walker_plot(samples, nwalkers, nsteps)
	fig = triangle.corner(samples, labels=[r'$a_{narrow}$', r'$mean_{narrow}$', r'$sigma_{narrow}$', r'$a_{broad}$', r'$mean_{broad}$', r'$sigma_{broad}$' ])
	fig.savefig('triangle_'+source[n]+'.pdf')
    	bf = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*N.percentile(samples, [16,50,84],axis=0)))
        N.save('best_fit_'+source[n]+'.npy', bf)
    	results = model([bf[0][0], bf[1][0], bf[2][0], bf[3][0], bf[4][0], bf[5][0]], wave)
    	print bf
    	P.figure()
    	P.plot(wave, flux, c='k', linewidth=2)
    	P.plot(wave, results, c='r', linestyle='dashed')
    	P.savefig('model_fit_'+source[n]+'.png')



