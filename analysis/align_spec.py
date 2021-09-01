#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 23:01:15 2020

@author: jlee
"""

import time
start_time = time.time()

import numpy as np
import glob, os
import copy
import sys
import init_cfg as ic
from astropy.io import fits
from scipy import interpolate


# ----- Loading the spatially aligned cubes ----- #
current_dir = os.getcwd()
dir_cb1 = ic.dir_cmb+"cube1/"
dir_cb2 = ic.dir_cmb+"cube2/"

working_dir = ic.dir_cmb
os.chdir(working_dir)

os.system("rm -rfv "+dir_cb2)
os.system("mkdir "+dir_cb2)


# ----- Making spectral combined images ----- #
wav_start, wav_end = ic.wav_range[0]+10., ic.wav_range[1]-10.

for i in np.arange(len(ic.cube_list)):
	hd0 = fits.getheader(dir_cb1+"al1_"+ic.cube_name[i]+"_3D.fits", ext=0)
	d_sci, h_sci = fits.getdata(dir_cb1+"al1_"+ic.cube_name[i]+"_3D.fits", ext=1, header=True)
	d_var, h_var = fits.getdata(dir_cb1+"al1_"+ic.cube_name[i]+"_3D.fits", ext=2, header=True)
	wav = np.linspace(start=h_sci['CRVAL3']+(1-h_sci['CRPIX3'])*h_sci['CD3_3'],
                      stop=h_sci['CRVAL3']+(h_sci['NAXIS3']-h_sci['CRPIX3'])*h_sci['CD3_3'],
                      num=h_sci['NAXIS3'], endpoint=True)

	d_scin = np.zeros((1+int((wav_end-wav_start)/ic.wav_intv), np.shape(d_sci)[1], np.shape(d_sci)[2]))
	d_varn = np.zeros((1+int((wav_end-wav_start)/ic.wav_intv), np.shape(d_sci)[1], np.shape(d_sci)[2]))
	wav_new = np.linspace(start=wav_start, stop=wav_end,
                          num=1+int((wav_end-wav_start)/ic.wav_intv), endpoint=True)

	for x in np.arange(d_sci.shape[2]):
		for y in np.arange(d_sci.shape[1]):
			func_sci = interpolate.interp1d(wav, d_sci[:,y,x], kind='linear')
			func_var = interpolate.interp1d(wav, d_var[:,y,x], kind='linear')
			d_scin[:,y,x] = func_sci(wav_new)
			d_varn[:,y,x] = func_var(wav_new)

	os.chdir(dir_cb2)

	nhd0 = fits.PrimaryHDU()
	nhd1 = fits.ImageHDU()
	nhd2 = fits.ImageHDU()

	nhd0.header = hd0

	nhd1.data = d_scin
	nhd1.header = h_sci
	nhd1.header['CRVAL3'] = wav_start
	nhd1.header['CD3_3'] = ic.wav_intv
	nhd1.header['CDELT3'] = ic.wav_intv	

	nhd2.data = d_varn
	nhd2.header = h_var
	nhd2.header['CRVAL3'] = wav_start
	nhd2.header['CD3_3'] = ic.wav_intv
	nhd2.header['CDELT3'] = ic.wav_intv

	cb_hdu = fits.HDUList([nhd0, nhd1, nhd2])
	cb_hdu.writeto("al2_"+ic.cube_name[i]+"_3D.fits", overwrite=True)

	print("Written: al2_"+ic.cube_name[i]+"_3D.fits")		

	os.chdir(working_dir)

os.chdir(current_dir)

print(f"Final wavelength range: {wav_start:.1f} - {wav_end:.1f} AA")


# Printing the running time
print('--- %.4f seconds ---' %(time.time()-start_time))
