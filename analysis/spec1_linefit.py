#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 15:01:41 2020

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import glob, os
from matplotlib import pyplot as plt
from astropy.io import fits
# from linefit import linefit
from linefit3 import linefit
# from linefit import linear


# ----- Basic parameters ----- #
redshift = 0.3424
dir_vbin = 'vorbin/'
dir_lines = 'lines3/'
os.system('rm -rfv '+dir_lines+'*')
os.system('mkdir '+dir_lines+'check/')


# ----- Loading Voronoi binned data ----- #
vb = np.load(dir_vbin+'vorbin_array.npz')
# wav, sci, var, cont
data_vbin = fits.getdata(dir_vbin+'vbin.fits').astype('int')
nvbin = np.unique(data_vbin).size-1
data_bfac_Ha = fits.getdata('bfac_bin_Ha.fits')
data_bfac_Hb = fits.getdata('bfac_bin_Hb.fits')


# ----- Line fitting ----- #

for num_line in [0, 3, 4, 6]:#[0, 1, 2, 3, 4]:

    # Loading linefit class
    if (num_line == 6):
        data_bfac = data_bfac_Hb
    else:
        data_bfac = data_bfac_Ha
    l = linefit(vb['wav'], vb['sci'], vb['var'], vb['cont'], num_line, redshift, dir_lines,
                broad_component=True, data_vbin=data_vbin, data_bfac=data_bfac)

    # 2D data initialization
    for ln in np.arange(l.nlines):
        os.system('mkdir '+dir_lines+l.line_names[ln])
        for col in ['mu', 'sigma', 'flux', 'vsig', 'snr', 'snrpix', 'rchisq',
                    'e_mu', 'e_vsig']:
            arr = np.zeros_like(data_vbin, dtype='float')
            exec(col+f"_{ln:d} = arr")
            
    # MCMC fitting of all the bins
    for ibin in np.arange(nvbin):
        df = l.solve(ibin, check=True, nwalkers=50,
                     ndiscard=1000, nsample=1000,
                     fluct0=1.0e-7, fluct1=1.0e-7, fluct2=1.0e-4, broad_component=True)
    
        if (ibin == 0):
            rchisq0 = np.array([])
            for ln in np.arange(l.nlines):
                rchisq0 = np.append(rchisq0, df['rchisq'].values[ln])

        theta = df['sigma'].values[0]
        for ln in np.arange(l.nlines):
            theta = np.append(theta, np.log(df['mu'].values[ln]))
            theta = np.append(theta, df['flux'].values[ln])

        # Solutions to 2D data
        nvbin_region = (data_vbin == ibin)
        N_area = np.sum(nvbin_region)

        prior_cnd = (np.isinf(l.log_prior(theta, ibin)) == False)
        sigma_sk_cnd = ((np.abs(df['g1_sigma'].values[0]) < 0.5) & \
                        (np.abs(df['k1_sigma'].values[0]) < 0.5))

        if (prior_cnd & sigma_sk_cnd):
            exec(f"sigma_{ln:d}[nvbin_region] = df['sigma'].values[ln]")
            for ln in np.arange(l.nlines):
                sk_cnd = ((np.abs(df['g1_mu'].values[ln]) < 0.5) & \
                          (np.abs(df['g1_flux'].values[ln]) < 0.5) & \
                          (np.abs(df['k1_mu'].values[ln]) < 0.5) & \
                          (np.abs(df['k1_flux'].values[ln]) < 0.5))
                if sk_cnd:
                    exec(f"mu_{ln:d}[nvbin_region] = df['mu'].values[ln]")
                    exec(f"sigma_{ln:d}[nvbin_region] = df['sigma'].values[ln]")
                    exec(f"flux_{ln:d}[nvbin_region] = df['flux'].values[ln] / N_area")
                    exec(f"vsig_{ln:d}[nvbin_region] = df['vsig'].values[ln]")
                    exec(f"snr_{ln:d}[nvbin_region] = df['snr'].values[ln]")
                    exec(f"snrpix_{ln:d}[nvbin_region] = df['snr'].values[ln] / np.sqrt(N_area)")
                    exec(f"rchisq_{ln:d}[nvbin_region] = df['rchisq'].values[ln] / rchisq0[ln]")
                    exec(f"e_mu_{ln:d}[nvbin_region] = df['e_mu'].values[ln]")
                    exec(f"e_vsig_{ln:d}[nvbin_region] = df['e_vsig'].values[ln]")

    # Saving the results
    for ln in np.arange(l.nlines):
        for col in ['mu', 'sigma', 'flux', 'vsig', 'snr', 'snrpix', 'rchisq',
                    'e_mu', 'e_vsig']:
            exec("arr = "+col+f"_{ln:d}")
            fits.writeto(dir_lines+l.line_names[ln]+'/'+col+'_2D.fits',
                         arr, overwrite=True)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
