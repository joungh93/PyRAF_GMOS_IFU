#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 11:26:25 2020

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import glob, os
import init_cfg as ic
from astropy.io import fits
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
import copy
from scipy.special import erf
from astropy.convolution import convolve
from astropy.convolution import Gaussian1DKernel
import tqdm


# ----- Directory ----- #
dir_lines = '/data/jlee/DATA/Gemini/Programs/GN-2019A-Q-215/analysis/lines/'
eline_name = ['NII6548', 'Halpha', 'NII6584']

for l in eline_name:
	os.system('rm -rfv '+dir_lines+l)
	os.system('mkdir '+dir_lines+l)
	exec('dir_'+l+' = '+"'"+dir_lines+l+"'")

current_dir = os.getcwd()


# ----- Function ----- #
def gauss_cdf_scale(x, mu, sigma, flux_scale):
    dx = x[1] - x[0]
    v1 = erf((x-mu+0.5*dx)/(np.sqrt(2.0)*sigma))
    v2 = erf((x-mu-0.5*dx)/(np.sqrt(2.0)*sigma))
    return flux_scale*(v1-v2)/(2.0*dx)

def multi3_gauss_cdf_scale(x, *pars):
	# g1 = gauss_cdf_scale(x, pars[0], pars[1], pars[2])
	# g2 = gauss_cdf_scale(x, pars[0]-14.75, pars[3], pars[4])
	# g3 = gauss_cdf_scale(x, pars[0]+20.66, pars[3], pars[5])

	g1 = gauss_cdf_scale(x, pars[0], pars[1], pars[2])
	g2 = gauss_cdf_scale(x, pars[3], pars[4], pars[5])
	g3 = gauss_cdf_scale(x, pars[6], pars[7], pars[8])

	return g1+g2+g3


# ----- Basic parameters ----- #
redshift = 0.353
dist_lum = 1875.5e+6  # pc
c = 2.99792e+5  # km/s
wav_NII_6548, wav_Ha, wav_NII_6584 = 6549.86, 6564.61, 6585.27


# ----- Reading the cube ----- #
fin_cb = 'sfcube_3D.fits'

hd0 = fits.getheader(fin_cb, ext=0)
d_sci2, h_sci = fits.getdata(fin_cb, ext=1, header=True)
d_var2, h_var = fits.getdata(fin_cb, ext=2, header=True)

wav = np.linspace(start=h_sci['CRVAL3'],
                  stop=h_sci['CRVAL3']+(h_sci['NAXIS3']-1)*h_sci['CD3_3'],
                  num=h_sci['NAXIS3'], endpoint=True)

wav_rest = wav/(1.0+redshift)
d_sci2 *= (1.0+redshift)
d_var2 *= (1.0+redshift)**2.0


# ----- Line fitting of all the pixels ----- #

# Variables
cont_2D = np.zeros((d_sci2.shape[1], d_sci2.shape[2]))
rms_2D = np.zeros((d_sci2.shape[1], d_sci2.shape[2])) 
rchisq_2D = np.zeros((d_sci2.shape[1], d_sci2.shape[2]))

param_name = ['lmu_2D' ,'lsig_2D', 'vsig_2D', 'flx_2D', 'snr_2D',
              'e_lmu_2D' ,'e_lsig_2D', 'e_vsig_2D', 'e_flx_2D']

for i in param_name:
	for j in eline_name:
		exec(i+'_'+j+' = np.zeros((d_sci2.shape[1], d_sci2.shape[2]))')

# Fitting
n_fit = 100

wav_fit = [6540.0, 6600.0]
spx_fit = [np.abs(wav_rest-wav_fit[0]).argmin(), np.abs(wav_rest-wav_fit[1]).argmin()]

spx_cont_left = [np.abs(wav_rest-6310.).argmin(), np.abs(wav_rest-6350.).argmin()]
spx_cont_right = [np.abs(wav_rest-6670.).argmin(), np.abs(wav_rest-6690.).argmin()]

x_bin = wav_rest[1] - wav_rest[0]
x_fit = wav_rest[spx_fit[0]:spx_fit[1]+1]

for x in tqdm.trange(d_sci2.shape[2]):
	for y in np.arange(d_sci2.shape[1]):

		cont_array = np.concatenate([d_sci2[spx_cont_left[0]:spx_cont_left[1], y, x],
			                         d_sci2[spx_cont_right[0]:spx_cont_right[1], y, x]], axis=0)
		cont = np.mean(cont_array)

		mpopt = []
		for i in np.arange(n_fit):
			y_dat = np.random.normal(d_sci2[spx_fit[0]:spx_fit[1]+1, y, x], np.sqrt(d_var2[spx_fit[0]:spx_fit[1]+1, y, x]))
			y_fit = y_dat-cont

			flx_scale0 = np.abs(np.sum(y_dat-cont)*x_bin)

			try:
				popt, pcov = curve_fit(multi3_gauss_cdf_scale, x_fit, y_fit,
					                   [wav_NII_6548, 3.0, 0.1*flx_scale0,
					                    wav_Ha, 3.0, 0.6*flx_scale0,
					                    wav_NII_6584, 3.0, 0.3*flx_scale0],
					                   bounds = ([x_fit[0], x_bin, 0.0,
					                   	          x_fit[0], x_bin, 0.0,
					                   	          x_fit[0], x_bin, 0.0],
					                   	         [x_fit[-1], wav_fit[1]-wav_fit[0], flx_scale0+0.1,
					                   	          x_fit[-1], wav_fit[1]-wav_fit[0], flx_scale0+0.1,
					                   	          x_fit[-1], wav_fit[1]-wav_fit[0], flx_scale0+0.1]))
				mpopt.append(popt)

			except RuntimeError:
				mpopt.append(np.zeros(9))

			mpar = np.mean(mpopt, axis=0)
			e_mpar = np.std(mpopt, axis=0)

		for l in np.arange(len(eline_name)):
			exec('lmu_2D_'+eline_name[l]+'[y, x] = mpar[{0:d}]'.format(l*3))
			exec('e_lmu_2D_'+eline_name[l]+'[y, x] = e_mpar[{0:d}]'.format(l*3))

			exec('lsig_2D_'+eline_name[l]+'[y, x] = mpar[{0:d}]'.format(l*3+1))
			exec('e_lsig_2D_'+eline_name[l]+'[y, x] = e_mpar[{0:d}]'.format(l*3+1))

			exec('flx_2D_'+eline_name[l]+'[y, x] = mpar[{0:d}]'.format(l*3+2))
			exec('e_flx_2D_'+eline_name[l]+'[y, x] = e_mpar[{0:d}]'.format(l*3+2))

			exec('vsig_2D_'+eline_name[l]+'[y, x] = c*mpar[{0:d}]/mpar[{1:d}]'.format(l*3+1, l*3))
			exec('e_vsig_2D_'+eline_name[l]+'[y, x] = vsig_2D_'+eline_name[l]+'[y, x]'+ \
			     ' * np.sqrt((e_mpar[{0:d}]/mpar[{0:d}])**2.0 + (e_mpar[{1:d}]/mpar[{1:d}])**2.0)'.format(l*3, l*3+1))

			exec('spx_lo = np.abs(wav_rest-(lmu_2D_'+eline_name[l]+'[y, x]-3.0*lsig_2D_'+eline_name[l]+'[y, x])).argmin()')
			exec('spx_hi = np.abs(wav_rest-(lmu_2D_'+eline_name[l]+'[y, x]+3.0*lsig_2D_'+eline_name[l]+'[y, x])).argmin()')

			exec('snr_2D_'+eline_name[l]+'[y, x] = flx_2D_'+eline_name[l]+'[y, x] / np.sqrt(np.sum(d_var2[spx_lo:spx_hi+1, y, x]))')

		cont_2D[y, x] = cont

		y_cal = cont + multi3_gauss_cdf_scale(x_fit, *mpar)
		y_obs = d_sci2[spx_fit[0]:spx_fit[1]+1, y, x]
		y_var = d_var2[spx_fit[0]:spx_fit[1]+1, y, x]

		rms_2D[y, x] = np.sqrt(np.sum((y_obs-y_cal)**2.0)/len(y_obs))
		rchisq_2D[y, x] = np.sum((y_obs-y_cal)**2.0 / y_var) / (len(y_obs)-9)


# ----- Saving the results in FITS files ----- #
for l in eline_name:
	exec('os.chdir(dir_'+l+')')

	quant_name = ['cont_2D', 'rms_2D', 'rchisq_2D',
	              'lmu_2D_'+l, 'lsig_2D_'+l, 'vsig_2D_'+l, 'flx_2D_'+l,
	              'snr_2D_'+l]

	for i in quant_name:
		hd0 = fits.PrimaryHDU()
		hd1 = fits.ImageHDU()
		exec('hd1.data = '+i)
		hd1.header['EXTNAME'] = 'SCI'

		if (i in ['lmu_2D_'+l ,'lsig_2D_'+l, 'vsig_2D_'+l, 'flx_2D_'+l]):
			hd2 = fits.ImageHDU()
			exec('hd2.data = e_'+i)
			hd2.header['EXTNAME'] = 'ERR'
			hdu = fits.HDUList([hd0, hd1, hd2])
		else:
			hdu = fits.HDUList([hd0, hd1])

		hdu.writeto(i.split('_'+l)[0]+'.fits', overwrite=True)


os.chdir(current_dir)


# Printing the running time
print('\n')
print('--- %s seconds ---' %(time.time()-start_time))