#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 16:53:20 2019

@author: jlee
"""


import numpy as np
import glob, os
from astropy.io import fits


# ----- Directories ----- #

# Basic directories
cpath = os.path.abspath(".")+"/"
dir_root = os.path.abspath("../")+"/"    # Same as 'dir_iraf'
dir_redux = dir_root+"redux/"
dir_wav = sorted(glob.glob(dir_redux+"w*"))

# Figure directory
dir_fig = cpath+"diagram/"
if (glob.glob(dir_fig) == []):
	os.system("mkdir "+dir_fig)

# Combine directory
dir_cmb = cpath+"combine/"
if (glob.glob(dir_cmb) == []):
	os.system("mkdir "+dir_cmb)


# ----- Basic configurations ----- #
centwave = [l.split("/")[-1][1:-1] for l in dir_wav]
cube_list, cube_name = [], []
for i in centwave:
	cube_list += sorted(glob.glob(dir_redux+"w"+i+"0/*/*_3D.fits"))
for i in np.arange(len(cube_list)):
	cube_name.append(cube_list[i].split('/')[-1].split('cstxeqxbrg')[-1].split('_3D.fits')[0])
cube_spa_off = []    # Cubes with spatial offset
cube_ref = ''    # Reference cube
pixel_scale = 0.1    # arcsec/pixel


# ----- Wavelength setting ----- #
redshift = 0.3424    # Redshift of galaxy
wav_range_res = np.array([6520.0, 6600.0])    # H alpha wavelength range (rest-frame)
check_x = [15, 55]    # [xmin, xmax] for check (including object and some bad regions)
check_y = [5, 45]    # [ymin, ymax] for check (including object and some bad regions)

# Reading wavelength range
wav_0, wav_1, dwav = [], [], []
for i in np.arange(len(cube_list)):
	h_sci = fits.getheader(cube_list[i], ext=1)
	wav = np.linspace(start=h_sci['CRVAL3']+(1-h_sci['CRPIX3'])*h_sci['CD3_3'],
                      stop=h_sci['CRVAL3']+(h_sci['NAXIS3']-h_sci['CRPIX3'])*h_sci['CD3_3'],
                      num=h_sci['NAXIS3'], endpoint=True)
	wav_0.append(wav[10])
	wav_1.append(wav[-11])
	dwav.append(h_sci['CD3_3'])
wav_0, wav_1, dwav = np.array(wav_0), np.array(wav_1), np.array(dwav)

# Total wavelength range
wav_range = [10 * (1 + wav_0.max() // 10),
             10 * (wav_1.min() // 10)]
nw_cut = int(round((wav_range[1]-wav_range[0])/np.mean(dwav))) + 1 
wav_intv = 1.0    # the resulting wavelength interval

combine_mode = 'median'    # 'median' / 'clippedmean'
