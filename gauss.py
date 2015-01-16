"""File to use emcee to fit a double gaussian with a narrow and broad peak to emission data from GANDALF."""

import numpy as N
import pylab as P
import emcee 
import triangle

def model(theta, wave):
	theta = a, u, s, a1, u1, s1
	return a * N.exp(- (((wave - u)**2) /(s**2))) + a1 * N.exp(- (((wave - u1)**2) /(s1**2)))

def lnlike(theta, wave, flux, fluxerr):
	pred_f = model(theta, wave)
	return -0.5*N.log(2*N.pi*fluxerr**2)-0.5*((flux-pred_f)**2/fluxerr**2)

def lnprior(theta):
	theta = a, u, s, a1, u1, s1
	if a >= 0 && a1 >= 0 && u>= 0 && u1 >= 0 && s >=0 && s1>=0:
		return 0.0
	else:
		return -N.inf

def lnprob(theta, wave, flux, fluxerr):
	lp = lnprior(theta)
    if not N.isfinite(lp):
        return -N.inf
    return lp + lnlike(theta, wave, flux, fluxerr)

ndim = 6
nwalkers = 100
nsteps = 1000
burnin = 1000
start = [22, 7800, 10, 10, 7800, 50]

#source = ['Q078017', 'Q507953', 'Q464967', 'Q292635', 'Q186225']
source = ['Q464967']

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


def walker_plot(samples, nwalkers, limit):
    """ Plotting function to visualise the steps of the walkers in each parameter dimension for smooth and disc theta values. 
        
        :samples:
        Array of shape (nsteps*nwalkers, 4) produced by the emcee EnsembleSampler in the sample function.
        
        :nwalkers:
        The number of walkers that step around the parameter space used to produce the samples by the sample function. Must be an even integer number larger than ndim.
        
        :limit:
        Integer value less than nsteps to plot the walker steps to. 
        
        RETURNS:
        :fig:
        The figure object
        """
    s = samples.reshape(nwalkers, -1, 4)
    s = s[:,:limit, :]
    fig = P.figure(figsize=(8,10))
    ax1 = P.subplot(4,1,1)
    ax2 = P.subplot(4,1,2)
    ax3 = P.subplot(4,1,3)
    ax4 = P.subplot(4,1,4)
    for n in range(len(s)):
        ax1.plot(s[n,:,0], 'k')
        ax2.plot(s[n,:,1], 'k')
        ax3.plot(s[n,:,2], 'k')
        ax4.plot(s[n,:,3], 'k')
    ax1.tick_params(axis='x', labelbottom='off')
    ax2.tick_params(axis='x', labelbottom='off')
    ax3.tick_params(axis='x', labelbottom='off')
    ax4.set_xlabel(r'step number')
    ax1.set_ylabel(r'$t_{smooth}$')
    ax2.set_ylabel(r'$\tau_{smooth}$')
    ax3.set_ylabel(r'$t_{disc}$')
    ax4.set_ylabel(r'$\tau_{disc}$')
    P.subplots_adjust(hspace=0.1)
    save_fig = 'walkers_steps_'+str(time.strftime('%H_%M_%d_%m_%y'))+'.pdf'
    fig.savefig(save_fig)
    return fig