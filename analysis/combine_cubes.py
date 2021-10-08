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
from astropy.io import fits
import tqdm


# ---- Loading the aligned cubes ----- #
current_dir = os.getcwd()
dir_cb2 = ic.dir_cmb+"cube2/"
working_dir = dir_cb2
os.chdir(working_dir)


# ----- Combining the spatially & spectrally aligned cubes ----- #
al2_cubes = sorted(glob.glob("al2_*_3D.fits"))
# p = np.genfromtxt(current_dir+"/count_max.txt", dtype=None, encoding=None, 
# 	              names=("cb","val"))
d_sci, h_sci = fits.getdata(al2_cubes[0], ext=1, header=True)
f_sci = np.zeros((d_sci.shape[0], d_sci.shape[1], d_sci.shape[2]))
f_var = np.zeros((d_sci.shape[0], d_sci.shape[1], d_sci.shape[2]))

for k in tqdm.trange(d_sci.shape[0]):
	ds = np.zeros((len(al2_cubes), d_sci.shape[1], d_sci.shape[2]))
	dv = np.zeros((len(al2_cubes), d_sci.shape[1], d_sci.shape[2]))
	if (len(ic.cube_spa_off) > 0):
		ds2 = np.zeros((len(al2_cubes)-len(ic.cube_spa_off), d_sci.shape[1], d_sci.shape[2]))

	j = 0
	for i in np.arange(len(al2_cubes)):
		hd0 = fits.getheader(al2_cubes[i], ext=0)
		d_sci, h_sci = fits.getdata(al2_cubes[i], ext=1, header=True)
		d_var, h_var = fits.getdata(al2_cubes[i], ext=2, header=True)
		
		# idx_cube = (p['cb'] == ic.cube_name[i])
		# count_max = p['val'][idx_cube]
		# d_sci[d_sci > count_max] = 0.0

		ds[i,:,:] = d_sci[k,:,:]
		dv[i,:,:] = d_var[k,:,:]

		if ((len(ic.cube_spa_off) > 0) & \
			(al2_cubes[i].split("al2_")[-1].split("_3D.fits")[0] not in ic.cube_spa_off)):
			ds2[j,:,:] = d_sci[k,:,:]
			j += 1

	if (ic.combine_mode == 'median'):
		cds = np.nanmedian(ds, axis=0)
	if (ic.combine_mode == 'clippedmean'):
		clipped_ds = sigma_clip(ds, sigma=3.0, axis=0).data
		cds = np.nanmean(clipped_ds, axis=0)
	cdv = np.sum(dv, axis=0) / ds.shape[0]**2.0

	##### Needed to be revised manually after visual check! #####
	if ((len(ic.cube_spa_off) > 0) & (ic.overwrite == True)):
		if (ic.combine_mode == 'median'):
			cds[:,61:] = np.median(ds2[:,:,61:], axis=0)
		if (ic.combine_mode == 'clippedmean'):
			clipped_ds2 = sigma_clip(ds2, sigma=3.0, axis=0).data
			cds[:,61:] = np.nanmean(clipped_ds2[:,:,61:], axis=0)
	##########

	f_sci[k,:,:] = cds
	f_var[k,:,:] = cdv


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
fcb_hdu.writeto("fcube_3D.fits", overwrite=True)


os.system("rm -rfv "+ic.dir_cmb+"split/")


# Printing the running time
print('--- %.4f seconds ---' %(time.time()-start_time))
