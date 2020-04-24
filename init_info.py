#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 16:53:20 2019

@author: jlee
"""


import numpy as np
import glob, os
from astropy.io import fits


current_dir = os.getcwd()
dir_raw = 'Raw/'


# ----- Collecting & sorting raw files ----- #
os.chdir(dir_raw)
rawfile = sorted(glob.glob('*.fits'))


# ----- Reading FITS headers ----- #
f = open(current_dir+'/'+'info.txt','w')
for i in rawfile:
	hdr = fits.getheader(i)
	objtype = hdr['OBSTYPE'].strip()
	objclass = hdr['OBSCLASS'].strip()
	centwave = str(hdr['CENTWAVE'])
	f.write(i.strip('.fits')+'   '+objtype+'   '+objclass+'   '+centwave+'\n')
f.close()

os.chdir(current_dir)


