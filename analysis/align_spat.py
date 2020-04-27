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
import init_cfg as ic


# ----- Importing IRAF from the root directory ----- #
current_dir = os.getcwd()
os.chdir(ic.dir_iraf)

from pyraf import iraf


# ----- Loading the processed cubes ----- #
from astropy.io import fits

dir_cmb = '/data/jlee/DATA/Gemini/Programs/GN-2019A-Q-215/analysis/combine/'
split_dir = dir_cmb+'split/'

working_dir = dir_cmb
os.chdir(working_dir)
iraf.chdir(working_dir)

os.system('rm -rfv '+split_dir)
os.system('mkdir '+split_dir)


# ----- Making the split 2D images for spatial alignment ----- #
wav_range = [5000, 9500]    # Total wavelength range
# check_Ha = [[2159,2170], [2168,2179]]

for i in np.arange(len(ic.cube_list)):
	fits_id = ic.cube_list[i].split('/')[-1].split('cstxeqxbrg')[-1].split('_3D.fits')[0]

	hd0 = fits.getheader(ic.cube_list[i], ext=0)
	d_sci, h_sci = fits.getdata(ic.cube_list[i], ext=1, header=True)
	d_var, h_var = fits.getdata(ic.cube_list[i], ext=2, header=True)
	# d_dq, h_dq = fits.getdata(ic.cube_list[i], ext=3, header=True)

	wav = np.linspace(start=h_sci['CRVAL3'], stop=h_sci['CRVAL3']+(h_sci['NAXIS3']-1)*h_sci['CD3_3'],
	                  num=h_sci['NAXIS3'], endpoint=True)

	# # Removing negative values
	# d_sci[d_sci < 0.0] = 0.0

	# # Removing bad pixel values
	# if (hd0['CENTWAVE'] == 700.0):
	# 	count_max = np.max(d_sci[check_Ha[0][0]:check_Ha[0][1], 2:46, 33:65])
	# if (hd0['CENTWAVE'] == 680.0):
	# 	count_max = np.max(d_sci[check_Ha[1][0]:check_Ha[1][1], 2:46, 33:65])

	# d_sci[d_sci > count_max] = 0.0

	# Wavelength cut
	wav_start, wav_end = wav_range[0], wav_range[1]
	nw_cut = int(round((wav_end-wav_start)/h_sci['CD3_3']))
	spx_start = np.abs(wav-wav_start).argmin()
	spx_end = spx_start + nw_cut

	d_sci_cut = d_sci[spx_start:spx_end,:,:]
	d_var_cut = d_var[spx_start:spx_end,:,:]
	# d_dq_cut = d_dq[spx_start:spx_end,:,:]

	os.chdir(split_dir)
	print('Working on '+fits_id+'\n')
	for j in np.arange(np.shape(d_sci_cut)[0]):
		fits.writeto(fits_id+'_2D_{0:04d}.fits'.format(j+1), d_sci_cut[j,:,:], h_sci, overwrite=True)
		fits.writeto(fits_id+'_VAR_{0:04d}.fits'.format(j+1), d_var_cut[j,:,:], h_var, overwrite=True)
	os.chdir(working_dir)


# Stop point #1
sys.exit('Split images generated')


# Shifting cubes with offset.txt (IRAF/imshift task)
os.chdir(split_dir)
iraf.unlearn('imshift')
for j in np.arange(np.shape(d_sci_cut)[0]):
	print('Shifting frame {0:04d}'.format(j+1))
	os.system('rm -rfv al1_N*_2D_{0:04d}.fits'.format(j+1))
	os.system('ls -1 N*_2D_{0:04d}.fits > split_{1:04d}.lis'.format(j+1, j+1))
	iraf.sleep(1.0)
	iraf.imshift(input='@split_{0:04d}.lis'.format(j+1), output='al1_//@split_{0:04d}.lis'.format(j+1),
		         shifts_file='../offset.txt', interp_type='linear', boundary_type='nearest')

	os.system('rm -rfv al1_N*_VAR_{0:04d}.fits'.format(j+1))
	os.system('ls -1 N*_VAR_{0:04d}.fits > variance_{1:04d}.lis'.format(j+1, j+1))
	iraf.sleep(1.0)
	iraf.imshift(input='@variance_{0:04d}.lis'.format(j+1), output='al1_//@variance_{0:04d}.lis'.format(j+1),
		         shifts_file='../offset.txt', interp_type='linear', boundary_type='nearest')

os.chdir(working_dir)


# Stop point #2
sys.exit('Shifted images generated.')



# Saving new cubes ([SCI, VAR])
os.system('rm -rfv '+working_dir+'cube1')
os.system('mkdir '+working_dir+'cube1')
for i in np.arange(len(ic.cube_list)):
	fits_id = ic.cube_list[i].split('/')[-1].split('cstxeqxbrg')[-1].split('_3D.fits')[0]

	hd0 = fits.getheader(ic.cube_list[i], ext=0)
	d_sci, h_sci = fits.getdata(ic.cube_list[i], ext=1, header=True)
	d_var, h_var = fits.getdata(ic.cube_list[i], ext=2, header=True)

	wav = np.linspace(start=h_sci['CRVAL3'], stop=h_sci['CRVAL3']+(h_sci['NAXIS3']-1)*h_sci['CD3_3'],
                      num=h_sci['NAXIS3'], endpoint=True)
	wav_start, wav_end = wav_range[0], wav_range[1]
	nw_cut = int(round((wav_end-wav_start)/h_sci['CD3_3']))
	spx_start = np.abs(wav-wav_start).argmin()
	spx_end = spx_start + nw_cut

	al1_d_sci = np.zeros((np.shape(d_sci_cut)[0], np.shape(d_sci_cut)[1], np.shape(d_sci_cut)[2]))
	al1_d_var = np.zeros((np.shape(d_sci_cut)[0], np.shape(d_sci_cut)[1], np.shape(d_sci_cut)[2]))

	for j in np.arange(np.shape(d_sci_cut)[0]):
		d_split = fits.getdata(split_dir+'al1_'+fits_id+'_2D_{0:04d}.fits'.format(j+1))
		al1_d_sci[j,:,:] = d_split
		d_var2 = fits.getdata(split_dir+'al1_'+fits_id+'_VAR_{0:04d}.fits'.format(j+1))
		al1_d_var[j,:,:] = d_var2

	os.chdir(working_dir+'cube1')

	nhd0 = fits.PrimaryHDU()
	nhd1 = fits.ImageHDU()
	nhd2 = fits.ImageHDU()

	nhd0.header = hd0

	nhd1.data = al1_d_sci
	nhd1.header = h_sci
	nhd1.header['CRVAL3'] = wav[spx_start]

	nhd2.data = al1_d_var
	nhd2.header = h_var
	nhd2.header['CRVAL3'] = wav[spx_start]

	cb_hdu = fits.HDUList([nhd0, nhd1, nhd2])
	cb_hdu.writeto('al1_'+fits_id+'_3D.fits', overwrite=True)

	print('Written : al1_'+fits_id+'_3D.fits')		

	os.chdir(working_dir)


os.chdir(current_dir)
iraf.chdir(current_dir)



# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))