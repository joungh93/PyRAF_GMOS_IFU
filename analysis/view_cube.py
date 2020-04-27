#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 16:53:20 2019

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
from matplotlib import pyplot as plt
import glob, os
from astropy.io import fits
import init_cfg as ic


# ----- Directories ----- #
dir_fig = '/data/jlee/DATA/Gemini/Programs/GN-2019A-Q-215/analysis/diagram/'
os.system('mkdir '+dir_fig+'cubes/')
dir_cmb = '/data/jlee/DATA/Gemini/Programs/GN-2019A-Q-215/analysis/combine/'


# ----- Wavelength cutout ranges ----- #
wav_range = [[6575, 6595],  # H beta
             [6760, 6790],  # [OIII] 5007
             [8520, 8540],  # [OI] 6300
             [8860, 8900],  # H alpha
             [8905, 8925],  # [NII] 6584
             [9080, 9130]]  # [SII] 6717,6731

check_Ha = [[2159,2170], [2168,2179]]
check_bad = [[2174,2176], [2200,2203]]


for i in np.arange(len(ic.cube_list)):
	hd0 = fits.getheader(ic.cube_list[i], ext=0)
	d_sci, h_sci = fits.getdata(ic.cube_list[i], ext=1, header=True)
	# d_var, h_var = fits.getdata(ic.cube_list[i], ext=2, header=True)

	wav = np.linspace(start=h_sci['CRVAL3'], stop=h_sci['CRVAL3']+(h_sci['NAXIS3']-1)*h_sci['CD3_3'],
	                  num=h_sci['NAXIS3'], endpoint=True)

	d_sci[d_sci < 0.0] = 0.0

	fits_id = ic.cube_list[i].split('/')[-1].split('cstxeqxbrg')[-1].split('_3D.fits')[0]

	# ----- Figure 1 : Histogram ----- #
	fig1 = plt.figure(i+1, figsize=(9,6))
	ax1 = plt.subplot(1,1,1)
	ax1.set_position([0.12,0.15,0.85,0.82])
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	ax1.set_xticks([1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0])
	ax1.set_xticklabels([r'$10^{-7}$', r'$10^{-6}$', r'$10^{-5}$', r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$'], fontsize=16.0)
	ax1.set_yticks([1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0])
	ax1.set_yticklabels([r'$10^{0}$', r'$10^{1}$', r'$10^{2}$', r'$10^{3}$', r'$10^{4}$', r'$10^{5}$'], fontsize=16.0)
	ax1.set_xlabel('Pixel values', fontsize=16.0)
	ax1.set_ylabel(r'$N$', fontsize=16.0)
	ax1.set_xlim([1.0e-7,2.0])
	ax1.set_ylim([9.0e-1, 5.0e+5])
	ax1.tick_params(width=2.0, length=12.0)
	plt.minorticks_on()
	ax1.tick_params(width=2.0, length=8.0, which='minor')
	for axis in ['top','bottom','left','right']:
	    ax1.spines[axis].set_linewidth(2.0)
	# -------------------------------- #

	# count, base = np.histogram(d_sci[:, 2:46, 33:65], bins=np.logspace(-7.0, 0.0, 70))
	p1 = ax1.hist(d_sci[:, 2:46, 33:65].ravel(), bins=np.logspace(-7.0, 0.0, 70), color='blue', alpha=0.75)
	
	if (hd0['CENTWAVE'] == 700.0):
		# count1, base = np.histogram(d_sci[check_Ha[0][0]:check_Ha[0][1], 2:46, 33:65], bins=np.logspace(-7.0, 0.0, 70))
		# count2, base = np.histogram(d_sci[check_bad[0][0]:check_bad[0][1], 2:46, 33:65], bins=np.logspace(-7.0, 0.0, 70))
		p2 = ax1.hist(d_sci[check_Ha[0][0]:check_Ha[0][1], 2:46, 33:65].ravel(), bins=np.logspace(-7.0, 0.0, 70), color='green', alpha=0.75)
		p3 = ax1.hist(d_sci[check_bad[0][0]:check_bad[0][1], 2:46, 33:65].ravel(), bins=np.logspace(-7.0, 0.0, 70), color='red', alpha=0.75)
		count_max = np.max(d_sci[check_Ha[0][0]:check_Ha[0][1], 2:46, 33:65])

	if (hd0['CENTWAVE'] == 680.0):
		# count1, base = np.histogram(d_sci[check_Ha[1][0]:check_Ha[1][1], 2:46, 33:65], bins=np.logspace(-7.0, 0.0, 70))
		# count2, base = np.histogram(d_sci[check_bad[1][0]:check_bad[1][1], 2:46, 33:65], bins=np.logspace(-7.0, 0.0, 70))
		p2 = ax1.hist(d_sci[check_Ha[1][0]:check_Ha[1][1], 2:46, 33:65].ravel(), bins=np.logspace(-7.0, 0.0, 70), color='green', alpha=0.75)
		p3 = ax1.hist(d_sci[check_bad[1][0]:check_bad[1][1], 2:46, 33:65].ravel(), bins=np.logspace(-7.0, 0.0, 70), color='red', alpha=0.75)
		count_max = np.max(d_sci[check_Ha[1][0]:check_Ha[1][1], 2:46, 33:65])

	ax1.plot([count_max, count_max], [9.0e-1, 1.0e+5], '--', color='k', linewidth=2.0, alpha=0.8)

	d_sci[d_sci > count_max] = 0.0

	plt.savefig(dir_fig+'cubes/'+'hist-'+fits_id+'.pdf')
	plt.close()


	d_col = np.zeros((np.shape(d_sci)[1], np.shape(d_sci)[2]))


	for j in np.arange(len(wav_range)):

		# ----- Figure 2 : wavelength cutout images ----- #
		fig1 = plt.figure(i+20+1, figsize=(10.2, 7.35))
		ax1 = plt.subplot(1,1,1)
		ax1.set_position([0.,0.,1.,1.])
		ax1.tick_params(labelleft=False)
		ax1.tick_params(labelbottom=False)
		ax1.tick_params(width=0.0, length=0.0)
		# ----------------------------------------------- #
		wav_start, wav_end = wav_range[j][0], wav_range[j][1]
		spx_start = np.abs(wav-wav_start).argmin()-1
		spx_end = np.abs(wav-wav_end).argmin()-1

		d_white = d_sci[spx_start:spx_end,:,:].sum(axis=0)
		v_low, v_high = np.percentile(d_white, [5.0, 95.0])

		ax1.imshow(d_white, cmap='gray', vmin=v_low, vmax=v_high, origin='lower')

		plt.savefig(dir_fig+'cubes/'+'sum-'+fits_id+'_'+str(wav_start)+'-'+str(wav_end)+'A.pdf')
		plt.close()

		d_col += d_white

		# Saving H alpha FITS images
		if (j == 3):
			fits.writeto(dir_cmb+'Ha_sum-'+fits_id+'.fits', d_white, h_sci, overwrite=True)


	# ----- Figure 3 : collapsed images ----- #
	fig1 = plt.figure(100+i, figsize=(10.2, 7.35))
	ax1 = plt.subplot(1,1,1)
	ax1.set_position([0.,0.,1.,1.])
	ax1.tick_params(labelleft=False)
	ax1.tick_params(labelbottom=False)
	ax1.tick_params(width=0.0, length=0.0)
	# --------------------------------------- #
	v_low, v_high = np.percentile(d_col, [5.0, 95.0])

	ax1.imshow(d_col, cmap='gray', vmin=v_low, vmax=v_high, origin='lower')

	plt.savefig(dir_fig+'cubes/'+'collpased-'+fits_id+'.pdf')
	plt.close()


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))