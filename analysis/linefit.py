#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 11:32:41 2020

@author: jlee
"""


import numpy as np
import copy
from astropy.convolution import convolve
from astropy.convolution import Gaussian1DKernel
from scipy.special import erf
from scipy.stats import sigmaclip
from scipy.optimize import minimize
import emcee
import pandas as pd
import warnings
from astropy.cosmology import FlatLambdaCDM
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec 
from scipy.stats import truncnorm
from scipy.stats import skew
from scipy.stats import kurtosis
from scipy.optimize import curve_fit


# ----- Function ----- #
def gauss_cdf_scale(x, mu, sigma, flux_scale):
    dx = x[1] - x[0]
    v1 = erf((x-mu+0.5*dx)/(np.sqrt(2.0)*sigma))
    v2 = erf((x-mu-0.5*dx)/(np.sqrt(2.0)*sigma))
    return flux_scale*(v1-v2)/(2.0*dx)


# ----- Class ----- #
class linefit:
    
    def __init__(self, wavelength, binned_spectrum, binned_variance, binned_continuum,
                 line_numbers, redshift, dir_lines,
                 broad_component=False, data_vbin=None, data_bfac=None):#, data_gaussian=None):

        '''
        wavelength :
        A wavelength data (1D array) on the observer-frame.
        The array should be the data type of (n_wavelength,).

        binned_spectrum :
        A spectrum data (2D array) after running voronoi 2D binning.
        The array should be the data type of (n_wavelength, n_bin).

        binned_variance :
        A variance data (2D array) after running voronoi 2D binning.
        The array should be the data type of (n_wavelength, n_bin).

        binned_continuum :
        A continuum data (2D array) after running voronoi 2D binning & continuum fitting.
        The array should be the data typde of (n_wavelength, n_bin).

        line_numbers : 
        A number of line declaration
        0 : [OII]3727/3729 line
        1 : H beta line
        2 : [OIII]4959/5007 line
        3 : H alpha + [NII]6548/6584 line
        4 : [SII]6717/6731 line
        5 : [OI]6300 line
        '''

        # Basic settings
        cosmo = FlatLambdaCDM(H0=70.0, Om0=0.3, Tcmb0=2.725)
        self.dir_lines = dir_lines
        self.redshift = redshift
        self.lumdist = cosmo.luminosity_distance(self.redshift).value * 1.0e+6  # pc
        self.c = 2.99792e+5  # km/s
        warnings.filterwarnings("ignore", category=RuntimeWarning) 


        # Reading the results of the integrated spectra
        fit_itg = np.genfromtxt('linefit_integrated.txt', dtype=None, encoding='ascii', comments='#',
                                names=('line','mu','e_mu','lsig','e_lsig','vsig','e_vsig',
                                       'R','e_R','flux','e_flux','rchisq'))
        self.fit_itg = fit_itg

        if broad_component:
            self.data_vbin = data_vbin
            self.data_bfac = data_bfac
            # self.data_gaussian = data_gaussian
            fit_itgb = np.genfromtxt('linefit_integrated_broad.txt', dtype=None, encoding='ascii', comments='#',
                                     names=('line','mu','e_mu','lsig','e_lsig','vsig','e_vsig',
                                            'R','e_R','flux','e_flux','rchisq','flxsum_scale'))
            self.fit_itgb = fit_itgb

        # Line declarations
        self.line_num = line_numbers

        if (line_numbers == 0):
            lname = ['OII3727']
            self.nlines = len(lname)
            self.line_names = lname
            self.line_wav = [fit_itg['mu'][fit_itg['line'] == l].item() for l in lname]
            self.wav_fit = [self.line_wav[0]-10.0, self.line_wav[-1]+10.0]

        if (line_numbers == 1):
            lname = ['Hbeta']
            self.nlines = len(lname)
            self.line_names = lname
            self.line_wav = [fit_itg['mu'][fit_itg['line'] == l].item() for l in lname]
            if broad_component:
                self.wav_fit = [self.line_wav[0]-50.0, self.line_wav[-1]+50.0]
            else:
                self.wav_fit = [self.line_wav[0]-10.0, self.line_wav[-1]+10.0]

        if (line_numbers == 2):
            lname = ['OIII4959', 'OIII5007']
            self.nlines = len(lname)
            self.line_names = lname
            self.line_wav = [fit_itg['mu'][fit_itg['line'] == l].item() for l in lname]
            if broad_component:
                self.wav_fit = [self.line_wav[0]-20.0, self.line_wav[-1]+30.0]
            else:
                self.wav_fit = [self.line_wav[0]-10.0, self.line_wav[-1]+10.0]

        if (line_numbers == 3):
            lname = ['NII6548', 'Halpha', 'NII6584']
            self.nlines = len(lname)
            self.line_names = lname
            self.line_wav = [fit_itg['mu'][fit_itg['line'] == l].item() for l in lname]
            if broad_component:
                self.wav_fit = [self.line_wav[0]-50.0, self.line_wav[-1]+40.0]
            else:
                self.wav_fit = [self.line_wav[0]-10.0, self.line_wav[-1]+10.0]

        if (line_numbers == 4):
            lname = ['SII6717', 'SII6731']
            self.nlines = len(lname)
            self.line_names = lname
            self.line_wav = [fit_itg['mu'][fit_itg['line'] == l].item() for l in lname]
            self.wav_fit = [self.line_wav[0]-10.0, self.line_wav[-1]+10.0]

        if (line_numbers == 5):
            lname = ['OI6300']
            self.nlines = len(lname)
            self.line_names = lname
            self.line_wav = [fit_itg['mu'][fit_itg['line'] == l].item() for l in lname]
            self.wav_fit = [self.line_wav[0]-10.0, self.line_wav[-1]+10.0]

        # Data
        self.wav_obs = wavelength
        self.wav_res = self.wav_obs / (1.0+self.redshift)
        self.spx_fit = [np.abs(self.wav_res-self.wav_fit[0]).argmin(),
                        np.abs(self.wav_res-self.wav_fit[1]).argmin()]
        self.nbin = binned_spectrum.shape[1]

        # Continuum subtraction
        data0 = binned_spectrum - binned_continuum
        vari0 = binned_variance
        cont0 = binned_continuum

        self.dat = data0 * (1.0+self.redshift)
        self.var = vari0 * (1.0+self.redshift)**2.0
        self.cont = cont0 * (1.0+self.redshift)

        # Reading the spectral resolution fitting results
        par, e_par = np.loadtxt('relation_wav_R.txt').T
        self.par = par
        self.e_par = e_par

        # Velocity sigma range
        Rmax = self.par[0] + self.par[1]*self.wav_obs[self.spx_fit[1]]
        e_Rmax = np.sqrt(self.e_par[0]**2.0 + (self.e_par[1]*self.wav_obs[self.spx_fit[1]])**2.0)
        vsig0 = self.c / (2.0*np.sqrt(2.0*np.log(2.0))*Rmax)
        e_vsig0 = self.c*e_Rmax / (2.0*np.sqrt(2.0*np.log(2.0))*Rmax*Rmax)
        self.vsig0 = vsig0
        self.vsig_llim = 0.5*self.vsig0
        if broad_component:
            self.vsig_ulim = 5.0*self.vsig0
        else:
            self.vsig_ulim = 2.0*self.vsig0

        # lsig0 = self.wav_res[self.spx_fit[1]] / (2.0*np.sqrt(2.0*np.log(2.0))*Rmax)
        # e_lsig0 = self.wav_res[self.spx_fit[1]]*e_Rmax / (2.0*np.sqrt(2.0*np.log(2.0))*Rmax*Rmax)
        # # lsig_llim = lsig0 - 3.0*e_lsig0
        # self.lsig0 = lsig0
        # self.lsig_llim = 0.5*self.lsig0
        # if broad_component:
        #     self.lsig_ulim = 5.0*self.lsig0
        # else:
        #     self.lsig_ulim = 2.0*self.lsig0


    def model_func(self, theta, x):
        dx = x[1]-x[0]
        val = 0.
        for i in np.arange(self.nlines):
            v1 = erf((x-theta[2*i+1]+0.5*dx)/(np.sqrt(2.0)*theta[0]))
            v2 = erf((x-theta[2*i+1]-0.5*dx)/(np.sqrt(2.0)*theta[0]))
            val += theta[2*i+2]*(v1-v2)/(2.0*dx)
        return val


    def log_likelihood(self, theta, x, y, yerr):
        mod = self.model_func(theta, x)
        sigma2 = yerr**2
        return -0.5 * np.sum((y-mod)**2 / sigma2 + np.log(2*np.pi*sigma2))


    def log_prior(self, theta, ibin):

        # Basic conditions
        icnd = 0
        sig_cnd = ((theta[0] > self.vsig_llim/self.c) & (theta[0] < self.vsig_ulim/self.c))
        icnd += 1*sig_cnd

        spec_sum = np.sum(np.abs(self.dat[self.spx_fit[0]:self.spx_fit[1]+1, ibin])* \
                          (self.wav_res[1]-self.wav_res[0]))

        for i in np.arange(self.nlines):
            mu_cnd = ((theta[2*i+1] > np.log(self.wav_fit[0])) & \
                      (theta[2*i+1] < np.log(self.wav_fit[1])))
            flx_cnd = ((theta[2*i+2] > 0.) & \
                       (theta[2*i+2] < 2.0*spec_sum))
            icnd += (1*mu_cnd + 1*flx_cnd)

        if (icnd == 2*self.nlines+1):
            return_value = 0.
        else:
            return_value = -np.inf

        # Specific conditions
        gauss_pdf = lambda X, M, S: np.exp(-0.5*((X-M)/S)**2.)/(S*np.sqrt(2.*np.pi))

        # Line 0, 1, 5: [OII]3727/3729, H beta, [OI]6300 (# of parameters = 3)
        if (self.line_num in [0, 1, 5]):
            vsig_init = self.fit_itg['vsig'][self.fit_itg['line'] == self.line_names[0]].item()
            e_vsig_init = self.fit_itg['e_vsig'][self.fit_itg['line'] == self.line_names[0]].item()
            fprior_sigma = np.log(gauss_pdf(theta[0], vsig_init/self.c, 10.*e_vsig_init/self.c))

            mu_init = self.fit_itg['mu'][self.fit_itg['line'] == self.line_names[0]].item()
            e_mu_init = self.fit_itg['e_mu'][self.fit_itg['line'] == self.line_names[0]].item()
            fprior_mu = np.log(gauss_pdf(theta[1], np.log(mu_init), 0.5/mu_init))

            fprior_flx = 0.

        # Line 2: [OIII]4959/5007 (# of parameters = 5)
        if (self.line_num == 2):
            vsig_init = self.fit_itg['vsig'][self.fit_itg['line'] == self.line_names[1]].item()
            e_vsig_init = self.fit_itg['e_vsig'][self.fit_itg['line'] == self.line_names[1]].item()
            fprior_sigma = np.log(gauss_pdf(theta[0], vsig_init/self.c, 10.*e_vsig_init/self.c))

            fprior_mu = 0.
            mu_init_arr = np.array([])
            for j in np.arange(self.nlines):
                mu_init = self.fit_itg['mu'][self.fit_itg['line'] == self.line_names[j]].item()
                e_mu_init = self.fit_itg['e_mu'][self.fit_itg['line'] == self.line_names[j]].item()
                fprior_mu += np.log(gauss_pdf(theta[2*j+1], np.log(mu_init), 0.5/mu_init))
                mu_init_arr = np.append(mu_init_arr, mu_init)
            fprior_mu += np.log(gauss_pdf(theta[3]-theta[1], np.log(mu_init_arr[1]/mu_init_arr[0]),
                                          0.01*np.log(mu_init_arr[1]/mu_init_arr[0])))

            flx2_cnd = (theta[4]/theta[2] > 1.)
            if flx2_cnd:
                fprior_flx = np.log(gauss_pdf(theta[4]/theta[2], 3.0, 0.1))  # prior function for flux ratio
            else:
                fprior_flx = -np.inf

        # Line 3: H alpha + [NII]6548/6584 (# of parameters = 7)
        if (self.line_num == 3):
            vsig_init = self.fit_itg['vsig'][self.fit_itg['line'] == self.line_names[1]].item()
            e_vsig_init = self.fit_itg['e_vsig'][self.fit_itg['line'] == self.line_names[1]].item()
            fprior_sigma = np.log(gauss_pdf(theta[0], vsig_init/self.c, 10.*e_vsig_init/self.c))

            fprior_mu = 0.
            mu_init_arr = np.array([])
            for j in np.arange(self.nlines):
                mu_init = self.fit_itg['mu'][self.fit_itg['line'] == self.line_names[j]].item()
                e_mu_init = self.fit_itg['e_mu'][self.fit_itg['line'] == self.line_names[j]].item()
                fprior_mu += np.log(gauss_pdf(theta[2*j+1], np.log(mu_init), 0.5/mu_init))
                mu_init_arr = np.append(mu_init_arr, mu_init)
            fprior_mu += np.log(gauss_pdf(theta[3]-theta[1], np.log(mu_init_arr[1]/mu_init_arr[0]),
                                          0.01*np.log(mu_init_arr[1]/mu_init_arr[0])))
            fprior_mu += np.log(gauss_pdf(theta[5]-theta[3], np.log(mu_init_arr[2]/mu_init_arr[1]),
                                          0.01*np.log(mu_init_arr[2]/mu_init_arr[1])))            

            flx2_cnd = ((theta[4]/theta[2] > 1.) & (theta[4]/theta[6] > 1.) & \
                        (theta[6]/theta[2] > 1.))
            if flx2_cnd:
                fprior_flx = np.log(gauss_pdf(theta[6]/theta[2], 3.0, 0.1))  # prior function for flux ratio
            else:
                fprior_flx = -np.inf

        # Line 4: [SII]6717/6731 (# of parameters = 5)
        if (self.line_num == 4):
            vsig_init = self.fit_itg['vsig'][self.fit_itg['line'] == self.line_names[0]].item()
            e_vsig_init = self.fit_itg['e_vsig'][self.fit_itg['line'] == self.line_names[0]].item()
            fprior_sigma = np.log(gauss_pdf(theta[0], vsig_init/self.c, 10.*e_vsig_init/self.c))

            fprior_mu = 0.
            mu_init_arr = np.array([])
            for j in np.arange(self.nlines):
                mu_init = self.fit_itg['mu'][self.fit_itg['line'] == self.line_names[j]].item()
                e_mu_init = self.fit_itg['e_mu'][self.fit_itg['line'] == self.line_names[j]].item()            
                fprior_mu += np.log(gauss_pdf(theta[2*j+1], np.log(mu_init), 0.5/mu_init))
                mu_init_arr = np.append(mu_init_arr, mu_init)
            fprior_mu += np.log(gauss_pdf(theta[3]-theta[1], np.log(mu_init_arr[1]/mu_init_arr[0]),
                                          0.01*np.log(mu_init_arr[1]/mu_init_arr[0])))

            fprior_flx = 0.

        return_value += (fprior_sigma + fprior_mu + fprior_flx)
    
        return return_value


    def solve(self, ibin, check=False, nwalkers=64, ndiscard=5000, nsample=5000,
              fluct0=1.0e-3, fluct1=1.0e-3, fluct2=1.0e-3, broad_component=False):

        ndim = 2*self.nlines+1          

        # Initial settings
        nll = lambda *args: -self.log_likelihood(*args)
        # for i in np.arange(self.nbin):
        Xfit = self.wav_res[self.spx_fit[0]:self.spx_fit[1]+1]
        Yfit = self.dat[self.spx_fit[0]:self.spx_fit[1]+1, ibin]
        e_Yfit = np.sqrt(self.var[self.spx_fit[0]:self.spx_fit[1]+1, ibin])
        spec_sum = np.sum(np.abs(Yfit)*(self.wav_res[1]-self.wav_res[0]))

        # Broad component subtraction
        if broad_component:
            broad_sum = np.zeros_like(Yfit)
            bfac = self.data_bfac[self.data_vbin == ibin][0]

            if (self.fit_itgb.size == 1):
                bline = self.fit_itgb['line'].item()
                if (self.line_num == int(bline[5])):
                    bpar = [self.fit_itgb['mu'].item(), self.fit_itgb['lsig'].item(),
                            bfac*self.fit_itgb['flux'].item()]
                    broad_sum += gauss_cdf_scale(Xfit, *bpar)
                else:
                    broad_sum += 0.

            else:
                bline = self.fit_itgb['line']
                bool_broad = pd.Series(bline).str.startswith(f"Broad{self.line_num:d}_").values
                nbroad = np.sum(bool_broad)
                if (nbroad == 0):
                    broad_sum += 0.
                else:
                    for b in np.arange(nbroad):
                        comp = (bline == f"Broad{self.line_num:d}_{b+1:d}")
                        bpar = [self.fit_itgb['mu'][comp].item(), self.fit_itgb['lsig'][comp].item(),
                                bfac*self.fit_itgb['flux'][comp].item()]
                        broad_sum += gauss_cdf_scale(Xfit, *bpar)

            Yfit = self.dat[self.spx_fit[0]:self.spx_fit[1]+1, ibin] - broad_sum

        Xfit2 = np.log(Xfit)
        Yfit2 = Xfit * Yfit
        e_Yfit2 = Xfit * e_Yfit

        # Finding the initial guess
        initial = np.zeros(ndim)
        initial[0] = self.vsig0 / self.c
        for j in np.arange(self.nlines):
            initial[2*j+1] = np.log(self.fit_itg['mu'][self.fit_itg['line'] == self.line_names[j]].item())
            initial[2*j+2] = spec_sum

        # Running MCMC
        log_posterior = lambda theta, x, y, yerr: self.log_prior(theta, ibin) + self.log_likelihood(theta, x, y, yerr) 

        pos = np.zeros((nwalkers, ndim))
        np.random.seed(0)
        pos[:,0] = truncnorm.rvs((self.vsig_llim/self.c - initial[0]) / fluct0,
                                 (self.vsig_ulim/self.c - initial[0]) / fluct0,
                                 loc=initial[0], scale=fluct0, size=nwalkers)

        for i in np.arange(1, ndim, 1):
            if (i % 2 == 1):
                pos[:,i] = truncnorm.rvs((np.log(self.wav_fit[0]) - initial[i]) / fluct1,
                                         (np.log(self.wav_fit[1]) - initial[i]) / fluct1,
                                         loc=initial[i], scale=fluct1, size=nwalkers)
            elif (i % 2 == 0):
                pos[:,i] = truncnorm.rvs((0. - initial[i]) / fluct2,
                                         (2.0*spec_sum - initial[i]) / fluct2,
                                         loc=initial[i], scale=fluct2, size=nwalkers)
        # print(pos)

        # pos = soln.x + fluct*np.random.randn(nwalkers, ndim)
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior,
                                        args=(Xfit2, Yfit2, e_Yfit2))
        sampler.run_mcmc(pos, ndiscard+nsample, progress=True)
        flat_samples = sampler.get_chain(discard=ndiscard, flat=True)
                
        popt, perr, labels, g1, k1 = [], [], [], [], []
        for j in np.arange(ndim):
            mcmc = np.percentile(flat_samples[:, j],
                                 [50-34.15, 50, 50+34.15])
            popt.append(mcmc[1])
            perr.append(0.5*(mcmc[2]-mcmc[0]))
            labels.append(f'par{j+1:d}')
            g1.append(skew(flat_samples[:, j]))
            k1.append(kurtosis(flat_samples[:, j]))
        popt = np.array(popt)
        perr = np.array(perr)
        g1 = np.array(g1)
        k1 = np.array(k1)

        popt2 = []
        popt2.append(popt[0]*np.median(Xfit))
        for j in 1+np.arange(ndim-1):
            if (j % 2 == 1):
                popt2.append(np.exp(popt[j]))
            else:
                popt2.append(popt[j])
        popt2 = np.array(popt2)


        if check:
            # Histogram plot
            fig = plt.figure(1, figsize=(10,10))
            gs = GridSpec(3, 3, left=0.05, bottom=0.05, right=0.975, top=0.975,
                          height_ratios=[1.]*3, width_ratios=[1.]*3,
                          hspace=0.15, wspace=0.15)
            for k in np.arange(ndim):
                ax = fig.add_subplot(gs[k // 3, k % 3])
                Y = flat_samples[:, k]
                ax.hist(Y, bins=20, histtype='step', linewidth=2.5)
                ax.tick_params(axis='both', labelsize=14.0)
                ax.tick_params(labelleft=False)
                ax.axvline(popt[k], color='k', linestyle='--', linewidth=2.5, alpha=0.7)
                ax.text(0.05, 0.95, f"g1 = {g1[k]:.2f}", fontsize=16.0, fontweight='bold',
                        color='red', ha='left', va='top', transform=ax.transAxes)
                ax.text(0.05, 0.85, f"k1 = {k1[k]:.2f}", fontsize=16.0, fontweight='bold',
                        color='green', ha='left', va='top', transform=ax.transAxes)
            plt.savefig(self.dir_lines+f"check/line{self.line_num:d}_bin{ibin:d}.png", dpi=300)
            plt.close()

            # Spectra plot
            fig = plt.figure(2, figsize=(12,9))
            ax = fig.add_subplot(111)
            ax.set_position([0.15,0.15,0.80,0.80])
            ax.tick_params(axis='both', labelsize=20.0)
            ax.set_xlim([Xfit[0]-25.0, Xfit[-1]+25.0])
            Ydat = self.dat[self.spx_fit[0]:self.spx_fit[1]+1, ibin]
            ax.set_ylim([np.min(Ydat)-1.0*np.abs(np.max(Ydat)), 1.25*np.abs(np.max(Ydat))])
            ax.set_xlabel(r"Rest-frame wavelength [${\rm \AA}$]", fontsize=20.0)
            ax.set_ylabel(r"Flux [${\rm 10^{-15}~erg~cm^{-2}~s^{-1}~\AA^{-1}}$]", fontsize=20.0)
            ax.plot(self.wav_res, self.dat[:, ibin], linewidth=3.0, alpha=0.7)
            ax.plot(self.wav_res, self.model_func(popt2, self.wav_res), linewidth=3.0, alpha=0.7)
            resi = -0.7*np.abs(np.max(Ydat))+self.dat[:, ibin]-self.model_func(popt2, self.wav_res)
            
            if broad_component:
                broad_sum_totwav = np.zeros_like(self.wav_res)
                if (self.fit_itgb.size == 1):
                    broad_totwav = gauss_cdf_scale(self.wav_res, self.fit_itgb['mu'].item(),
                                                   self.fit_itgb['lsig'].item(), bfac*self.fit_itgb['flux'].item())
                    ax.plot(self.wav_res, broad_totwav,
                            linewidth=2.5, linestyle='--', color='red', alpha=0.6)
                    broad_sum_totwav += broad_totwav
                else:
                    for b in np.arange(self.fit_itgb.size):
                        broad_totwav = gauss_cdf_scale(self.wav_res, self.fit_itgb['mu'][b],
                                                       self.fit_itgb['lsig'][b], self.fit_itgb['flux'][b]*bfac)
                    ax.plot(self.wav_res, broad_totwav,
                            linewidth=2.5, linestyle='--', color='red', alpha=0.6)
                    broad_sum_totwav += broad_totwav    
                resi -= broad_sum_totwav
            
            ax.plot(self.wav_res, resi, linewidth=2.5, color='green', alpha=0.6)
            plt.savefig(self.dir_lines+f"check/line{self.line_num:d}_fit_bin{ibin:d}.png", dpi=300)
            plt.close()


        vsig = self.c*popt[0]
        e_vsig = self.c*perr[0]
        snr = popt[2::2] / perr[2::2]

        Ymod = self.model_func(popt, np.log(self.wav_res))[self.spx_fit[0]:self.spx_fit[1]+1]
        rchisq = []
        for j in np.arange(self.nlines):
            spx_line = [np.abs(Xfit2-(popt[1::2][j]-3*popt[0])).argmin(),
                        np.abs(Xfit2-(popt[1::2][j]+3*popt[0])).argmin()]
            chisq = ((Yfit2-Ymod)/e_Yfit2)**2.
            dof = len(Yfit2[spx_line[0]:spx_line[1]+1])-3
            rchisq.append(np.sum(chisq[spx_line[0]:spx_line[1]+1]) / dof)
        # rchisq = np.sum(((Yfit-Ymod)/e_Yfit)**2.) / (len(Yfit)-ndim)

        if broad_component:
            broad = bfac
        else:
            broad = 0.0

        df = pd.DataFrame(data = {'line': self.line_names,
                                  'mu': np.exp(popt[1::2]),#popt[0::3],
                                  'e_mu': np.exp(popt[1::2])*perr[1::2],#perr[0::3],
                                  'g1_mu': g1[1::2],
                                  'k1_mu': k1[1::2],
                                  'sigma': [popt[0]]*len(self.line_names),#popt[1::3],
                                  'e_sigma': [perr[0]]*len(self.line_names),#perr[1::3],
                                  'g1_sigma': [g1[0]]*len(self.line_names),
                                  'k1_sigma': [k1[0]]*len(self.line_names),
                                  'flux': popt[2::2],
                                  'e_flux': perr[2::2],
                                  'g1_flux': g1[2::2],
                                  'k1_flux': k1[2::2],
                                  'vsig': vsig,
                                  'e_vsig': e_vsig,
                                  'snr': snr,
                                  'rchisq': rchisq,
                                  'broad': broad})

        return df


if (__name__ == '__main__'):

    import numpy as np
    import glob, os
    from matplotlib import pyplot as plt
    from astropy.io import fits

    # ----- Basic parameters ----- #
    redshift = 0.3527
    dir_vbin = 'vorbin/'
    dir_lines = 'lines3/'
    if (glob.glob(dir_lines) == []):
        os.system("mkdir "+dir_lines)
        if (glob.glob(dir_lines+"check/") == []):
            os.system("mkdir "+dir_lines+"check/")


    # ----- Loading Voronoi binned data ----- #
    vb = np.load(dir_vbin+'vorbin_array.npz')
    # wav, sci, var
    data_vbin = fits.getdata(dir_vbin+'vbin.fits').astype('int')
    nvbin = np.unique(data_vbin).size-1

    l1 = linefit(vb['wav'], vb['sci'], vb['var'], vb['cont'], 1, redshift, dir_lines)
    l2 = linefit(vb['wav'], vb['sci'], vb['var'], vb['cont'], 2, redshift, dir_lines)
    l3 = linefit(vb['wav'], vb['sci'], vb['var'], vb['cont'], 3, redshift, dir_lines)
    l4 = linefit(vb['wav'], vb['sci'], vb['var'], vb['cont'], 4, redshift, dir_lines)

    test_ibin = [0, 1, 2, 3, 6, 22, 50, 53]

    for ibin in test_ibin:

        # df1 = l1.solve(ibin, check=True, nwalkers=50,
        #                ndiscard=1000, nsample=1000,
        #                fluct0=1.0e-6, fluct1=1.0e-7, fluct2=1.0e-4)
        # theta1 = df1['sigma'].values[0]
        # for ln in np.arange(l1.nlines):
        #     theta1 = np.append(theta1, np.log(df1['mu'].values[ln]))
        #     theta1 = np.append(theta1, df1['flux'].values[ln])
        # print(l1.log_prior(theta1, ibin))

        # df2 = l2.solve(ibin, check=True, nwalkers=50,
        #                ndiscard=1000, nsample=1000,
        #                fluct0=1.0e-6, fluct1=1.0e-7, fluct2=1.0e-4)
        # theta2 = df2['sigma'].values[0]
        # for ln in np.arange(l2.nlines):
        #     theta2 = np.append(theta2, np.log(df2['mu'].values[ln]))
        #     theta2 = np.append(theta2, df2['flux'].values[ln])
        # print(l2.log_prior(theta2, ibin))

        df3 = l3.solve(ibin, check=True, nwalkers=50,
                       ndiscard=1000, nsample=1000,
                       fluct0=1.0e-6, fluct1=1.0e-7, fluct2=1.0e-4)
        theta3 = df3['sigma'].values[0]
        for ln in np.arange(l3.nlines):
            theta3 = np.append(theta3, np.log(df3['mu'].values[ln]))
            theta3 = np.append(theta3, df3['flux'].values[ln])
        print(l3.log_prior(theta3, ibin))

        # df4 = l4.solve(ibin, check=True, nwalkers=50,
        #                ndiscard=1000, nsample=1000,
        #                fluct0=1.0e-6, fluct1=1.0e-7, fluct2=1.0e-4)
        # theta4 = df4['sigma'].values[0]
        # for ln in np.arange(l4.nlines):
        #     theta4 = np.append(theta4, np.log(df4['mu'].values[ln]))
        #     theta4 = np.append(theta4, df4['flux'].values[ln])
        # print(l4.log_prior(theta4, ibin))
