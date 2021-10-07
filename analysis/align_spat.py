#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 15:18:20 2020
@author: jlee
"""

import time
start_time = time.time()

import numpy as np
import glob, os
import copy
import sys
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from scipy import ndimage
from tqdm import trange
import init_cfg as ic
import warnings
warnings.filterwarnings("error")


# ----- Directories ----- #
current_dir = os.getcwd()
working_dir = ic.dir_cmb
split_dir = ic.dir_cmb+"split/"
os.chdir(working_dir)
os.system('rm -rfv '+split_dir)
os.system('mkdir '+split_dir)


# ----- Making the split 2D images for spatial alignment ----- #
kernel = Gaussian2DKernel(2)
p = np.genfromtxt(current_dir+"/count_max.txt", dtype=None, encoding=None, 
	              names=("cb","val"))
count_max = p['val'].max()

# Making the split images
for i in np.arange(len(ic.cube_list)):
	hd0 = fits.getheader(ic.cube_list[i], ext=0)
	d_sci, h_sci = fits.getdata(ic.cube_list[i], ext=1, header=True)
	d_var, h_var = fits.getdata(ic.cube_list[i], ext=2, header=True)

	wav = np.linspace(start=h_sci['CRVAL3']+(1-h_sci['CRPIX3'])*h_sci['CD3_3'],
                      stop=h_sci['CRVAL3']+(h_sci['NAXIS3']-h_sci['CRPIX3'])*h_sci['CD3_3'],
                      num=h_sci['NAXIS3'], endpoint=True)

	# Wavelength cut
	spx_start = np.abs(wav-ic.wav_range[0]).argmin()
	spx_end = spx_start + ic.nw_cut

	d_sci_cut = d_sci[spx_start:spx_end,:,:]
	d_var_cut = d_var[spx_start:spx_end,:,:]

	assert (p['cb'][i] == ic.cube_name[i])

	os.chdir(split_dir)
	print('Working on '+ic.cube_name[i])
	for j in trange(np.shape(d_sci_cut)[0]):
		d_sci_cut2 = d_sci_cut[j,:,:]
		d_var_cut2 = d_var_cut[j,:,:]

		d_sci_cut2[d_sci_cut2 > count_max] = np.nan
		d_sci_cut2[d_sci_cut2 < -count_max] = 0.0

		try:
			conv = convolve(d_sci_cut2, kernel)
			d_sci_cut2[np.isnan(d_sci_cut2)] = conv[np.isnan(d_sci_cut2)]
		except:
			try:
				conv2 = convolve(d_sci_cut2, Gaussian2DKernel(4))
				d_sci_cut2[np.isnan(d_sci_cut2)] = conv2[np.isnan(d_sci_cut2)]
			except:
				d_sci_cut2[np.isnan(d_sci_cut2)] = 0.5*(count_max + np.nanmedian(d_sci_cut2))

		fits.writeto(ic.cube_name[i]+f"_SCI_{j+1:04d}.fits", d_sci_cut2, h_sci, overwrite=True)
		fits.writeto(ic.cube_name[i]+f"_VAR_{j+1:04d}.fits", d_var_cut2, h_var, overwrite=True)
	os.chdir(working_dir)
print('----- Split images were created. -----\n')


# ----- Shifting cubes with offset.txt ----- #
os.chdir(split_dir)
offset_X, offset_Y = np.loadtxt(ic.dir_cmb+"offset.txt").T
for j in trange(d_sci_cut.shape[0]):
	for k in np.arange(len(ic.cube_name)):
		ds, hs = fits.getdata(ic.cube_name[k]+f"_SCI_{j+1:04d}.fits", header=True)
		ds_shifted = ndimage.shift(ds, shift=(-offset_Y[k], -offset_X[k]), mode='nearest')
		fits.writeto("al1_"+ic.cube_name[k]+f"_SCI_{j+1:04d}.fits", ds_shifted, hs, overwrite=True)

		dv, hv = fits.getdata(ic.cube_name[k]+f"_VAR_{j+1:04d}.fits", header=True)
		dv_shifted = ndimage.shift(dv, shift=(-offset_Y[k], -offset_X[k]), mode='nearest')
		fits.writeto("al1_"+ic.cube_name[k]+f"_VAR_{j+1:04d}.fits", dv_shifted, hv, overwrite=True)

os.chdir(working_dir)	
print('----- Shift images were created. -----\n')


# ----- Saving new cubes ([SCI, VAR]) ----- #
os.system("rm -rfv "+working_dir+"cube1")
os.system("mkdir "+working_dir+"cube1")

for i in np.arange(len(ic.cube_list)):
	hd0 = fits.getheader(ic.cube_list[i], ext=0)
	d_sci, h_sci = fits.getdata(ic.cube_list[i], ext=1, header=True)
	d_var, h_var = fits.getdata(ic.cube_list[i], ext=2, header=True)

	# assert (p['cb'][i] == ic.cube_name[i])
	# count_max = p['val'].max()
	# d_sci[d_sci > count_max] = 0.0
	
	wav = np.linspace(start=h_sci['CRVAL3']+(1-h_sci['CRPIX3'])*h_sci['CD3_3'],
                      stop=h_sci['CRVAL3']+(h_sci['NAXIS3']-h_sci['CRPIX3'])*h_sci['CD3_3'],
                      num=h_sci['NAXIS3'], endpoint=True)

	spx_start = np.abs(wav-ic.wav_range[0]).argmin()
	spx_end = spx_start + ic.nw_cut

	d_sci_cut = d_sci[spx_start:spx_end,:,:]
	d_var_cut = d_var[spx_start:spx_end,:,:]

	al1_d_sci = np.zeros((np.shape(d_sci_cut)[0], np.shape(d_sci_cut)[1], np.shape(d_sci_cut)[2]))
	al1_d_var = np.zeros((np.shape(d_sci_cut)[0], np.shape(d_sci_cut)[1], np.shape(d_sci_cut)[2]))

	for j in np.arange(np.shape(d_sci_cut)[0]):
		d_sci2 = fits.getdata(split_dir+"al1_"+ic.cube_name[i]+f"_SCI_{j+1:04d}.fits")
		al1_d_sci[j,:,:] = d_sci2
		d_var2 = fits.getdata(split_dir+"al1_"+ic.cube_name[i]+f"_VAR_{j+1:04d}.fits")
		al1_d_var[j,:,:] = d_var2

	os.chdir(working_dir+"cube1")

	nhd0 = fits.PrimaryHDU()
	nhd1 = fits.ImageHDU()
	nhd2 = fits.ImageHDU()

	nhd0.header = hd0

	nhd1.data = al1_d_sci
	nhd1.header = h_sci
	nhd1.header['CRPIX3'] = 1.
	nhd1.header['CRVAL3'] = wav[spx_start]

	nhd2.data = al1_d_var
	nhd2.header = h_var
	nhd2.header['CRPIX3'] = 1.
	nhd2.header['CRVAL3'] = wav[spx_start]

	cb_hdu = fits.HDUList([nhd0, nhd1, nhd2])
	cb_hdu.writeto("al1_"+ic.cube_name[i]+"_3D.fits", overwrite=True)

	print("Written: al1_"+ic.cube_name[i]+"_3D.fits")		

	os.chdir(working_dir)

os.chdir(current_dir)


# Printing the running time
print('--- %.4f seconds ---' %(time.time()-start_time))
