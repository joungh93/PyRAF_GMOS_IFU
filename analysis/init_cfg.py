#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 16:53:20 2019

@author: jlee
"""


import numpy as np
import glob, os


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


# ----- Wavelength setting ----- #
redshift = 0.3424    # Redshift of galaxy
wav_range = np.array([6520.0, 6600.0])    # H alpha wavelength range (rest-frame)
check_x = [15, 55]    # [xmin, xmax] for check
check_y = [5, 45]    # [ymin, ymax] for check

