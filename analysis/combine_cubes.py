#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 18:10:15 2020

@author: jlee
"""

import time
start_time = time.time()

import numpy as np
import glob, os
import init_cfg as ic
from astropy.stats import sigma_clip
import tqdm


# ---- Loading the aligned cubes ----- #
current_dir = os.getcwd()

dir_cmb = ic.dir_iraf+'analysis/combine/'
dir_cb2 = dir_cmb+'cube2/'

working_dir = dir_cb2
os.chdir(working_dir)


# ----- Combining the spatially & spectrally aligned cubes ----- #
from astropy.io import fits

al2_cubes = sorted(glob.glob('al2_*_3D.fits'))
d_sci, h_sci = fits.getdata(al2_cubes[0], ext=1, header=True)

f_sci = np.zeros((d_sci.shape[0], d_sci.shape[1], d_sci.shape[2]))
f_var = np.zeros((d_sci.shape[0], d_sci.shape[1], d_sci.shape[2]))

for k in tqdm.trange(d_sci.shape[0]):
	d1s = np.zeros((len(al2_cubes)-len(ic.cube_spa_off), d_sci.shape[1], d_sci.shape[2]))
	d1v = np.zeros((len(al2_cubes)-len(ic.cube_spa_off), d_sci.shape[1], d_sci.shape[2]))
	d2s = np.zeros((len(ic.cube_spa_off), d_sci.shape[1], d_sci.shape[2]))
	d2v = np.zeros((len(ic.cube_spa_off), d_sci.shape[1], d_sci.shape[2]))

	num_cb_noff = 0
	num_cb_off = 0

	for i in np.arange(len(al2_cubes)):
		cb = al2_cubes[i].split('_')[1]

		hd0 = fits.getheader(al2_cubes[i], ext=0)
		d_sci, h_sci = fits.getdata(al2_cubes[i], ext=1, header=True)
		d_var, h_var = fits.getdata(al2_cubes[i], ext=2, header=True)

		if (cb not in ic.cube_spa_off):
			d1s[num_cb_noff,:,:] = d_sci[k,:,:]
			d1v[num_cb_noff,:,:] = d_var[k,:,:]
			num_cb_noff += 1
		else:
			d2s[num_cb_off,:,:] = d_sci[k,:,:]
			d2v[num_cb_off,:,:] = d_var[k,:,:]
			num_cb_off += 1

	clipped_d1s = sigma_clip(d1s, sigma=3.0, axis=0).data
	cd1_sci = np.nanmean(clipped_d1s, axis=0)
	cd1_var = np.sum(d1v, axis=0) / num_cb_noff**2.0

	clipped_d2s = sigma_clip(d2s, sigma=3.0, axis=0).data
	cd2_sci = np.nanmean(clipped_d2s, axis=0)
	cd2_var = np.sum(d2v, axis=0) / num_cb_off**2.0

	f_sci[k,:,:] = 0.5*(cd1_sci+cd2_sci)
	f_sci[k,:,62-1:] = cd1_sci[:,62-1:]


	f_var[k,:,:] = 0.5*0.5*(cd1_var+cd2_var)
	f_var[k,:,62-1:] = cd1_var[:,62-1:]
	# f_var[k,:,:] = 0.5*(med1_var+med2_var)
	# f_var[k,:,62-1:] = med1_var[:,62-1:]

	# print("Frame {0:04d}/{1:04d} of the wavelength axis is created.".format(k+1, d_sci.shape[0]))


# ----- Creating the final cube for analysis ----- #
os.chdir(current_dir)

fhd0 = fits.PrimaryHDU()
fhd1 = fits.ImageHDU()
fhd2 = fits.ImageHDU()

fhd0.header = hd0

fhd1.data = f_sci
fhd1.header = h_sci

fhd2.data = f_var
fhd2.header = h_var

fcb_hdu = fits.HDUList([fhd0, fhd1, fhd2])
fcb_hdu.writeto('fcube_3D.fits', overwrite=True)


os.system('rm -rfv '+dir_cmb+'split/')


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
