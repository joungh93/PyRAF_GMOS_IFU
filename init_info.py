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
dir_raw = 'raw/'


# ----- Collecting & sorting raw files ----- #
os.chdir(dir_raw)
rawfile = sorted(glob.glob('*.fits'))


# ----- Reading FITS headers ----- #
f = open(current_dir+'/'+'info.txt','w')
f.write("# FILENAME  OBJTYPE  OBSCLASS  CENTWAVE  DATALAB  EXPTIME  MASKNAME  GRATING  AIRMASS  MJD-OBS\n")
for i in rawfile:
	h0 = fits.getheader(i, ext=0)
	h1 = fits.getheader(i, ext=1)
	objtype = h0['OBSTYPE'].strip()
	objclass = h0['OBSCLASS'].strip()
	centwave = str(h0['CENTWAVE'])
	datalabel = h0['DATALAB'].strip()
	exptime = f"{h0['EXPTIME']:.1f}"
	mask = h0['MASKNAME'].strip()
	grating = h0['GRATING'].strip()
	airmass = f"{h0['AIRMASS']:.4f}"
	mjd = f"{h1['MJD-OBS']:f}"

	f.write(i.strip('.fits')+'  ')
	f.write(objtype+'  ')
	f.write(objclass+'  ')
	f.write(centwave+'  ')
	f.write(datalabel+'  ')
	f.write(exptime+'  ')
	f.write(mask+'  ')
	f.write(grating+'  ')
	f.write(airmass+'  ')
	f.write(mjd+'\n')
f.close()

os.chdir(current_dir)
