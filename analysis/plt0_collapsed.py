#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 5 15:52:50 2020

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
from matplotlib import pyplot as plt
import glob, os
from astropy.io import fits
import init_cfg as ic


dir_fig = '/data/jlee/DATA/Gemini/Programs/GN-2019A-Q-215/analysis/diagram/'


# ----- Reading the final cube ----- #
dir_fig = ic.cpath+"diagram/"
fin_cb = "bfcube_3D.fits"
hd0 = fits.getheader(fin_cb, ext=0)
d_sci, h_sci = fits.getdata(fin_cb, ext=1, header=True)
d_var, h_var = fits.getdata(fin_cb, ext=2, header=True)
wav = np.linspace(start=h_sci['CRVAL3'],
                  stop=h_sci['CRVAL3']+(h_sci['NAXIS3']-1)*h_sci['CD3_3'],
                  num=h_sci['NAXIS3'], endpoint=True)


# ----- Making the collapsed images ----- #
# d_sci[d_sci < 0.0] = 0.0

### Figure 1 : white images
fig1 = plt.figure(1, figsize=(10.2, 7.35))
ax1 = plt.subplot(1,1,1)
ax1.set_position([0.,0.,1.,1.])
ax1.tick_params(labelleft=False)
ax1.tick_params(labelbottom=False)
ax1.tick_params(width=0.0, length=0.0)
# ----------------------- #
d_white = d_sci.sum(axis=0)
v_low, v_high = np.percentile(d_white, [5.0, 95.0])

ax1.imshow(d_white, cmap='gray', vmin=v_low, vmax=v_high, origin='lower')

plt.savefig(dir_fig+'continuum.pdf')
plt.close()


# ----- Making the wavelength cutout images ----- #

# Wavelength cutout ranges
wav_range = [[5035, 5055],  # [OII]3727,3729
             [5865, 5885],  # H gamma
             [6570, 6590],  # H beta
             [6700, 6720],  # [OIII]4959
             [6765, 6785],  # [OIII]5007
             [8520, 8535],  # [OI]6300
             [8870, 8895],  # H alpha
             [8900, 8917],  # [NII]6584
             [9080, 9120]]  # [SII]6717,6731
emissions = ['OII3727+29', 'Hgamma','Hbeta','OIII4959','OIII5007',
             'OI6300','Halpha','NII6584','SII6717+31']

d_col = np.zeros((np.shape(d_sci)[1], np.shape(d_sci)[2]))

for j in np.arange(len(emissions)):

	### Figure 2 : wavelength cutout images
	fig1 = plt.figure(j+10+1, figsize=(10.2, 7.35))
	ax1 = plt.subplot(1,1,1)
	ax1.set_position([0.,0.,1.,1.])
	ax1.tick_params(labelleft=False)
	ax1.tick_params(labelbottom=False)
	ax1.tick_params(width=0.0, length=0.0)
	# ----------------------------------- #
	wav_start, wav_end = wav_range[j][0], wav_range[j][1]
	spx_start = np.abs(wav-wav_start).argmin()
	spx_end = np.abs(wav-wav_end).argmin()

	d_lines = d_sci[spx_start:spx_end,:,:].sum(axis=0)
	v_low, v_high = np.percentile(d_lines, [5.0, 95.0])

	ax1.imshow(d_lines, cmap='gray', vmin=v_low, vmax=v_high, origin='lower')

	plt.savefig(dir_fig+'wavcut_'+emissions[j]+'.pdf')
	plt.close()

	d_col += d_lines


### Figure 3 : emission line images
fig1 = plt.figure(50, figsize=(10.2, 7.35))
ax1 = plt.subplot(1,1,1)
ax1.set_position([0.,0.,1.,1.])
ax1.tick_params(labelleft=False)
ax1.tick_params(labelbottom=False)
ax1.tick_params(width=0.0, length=0.0)
# ------------------------------- #
v_low, v_high = np.percentile(d_col, [5.0, 95.0])

ax1.imshow(d_col, cmap='gray', vmin=v_low, vmax=v_high, origin='lower')

plt.savefig(dir_fig+'emissions.pdf')
plt.close()


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
