#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 3 01:18:27 2020

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import glob, os
import g0_init_cfg as ic


# ----- Importing IRAF from the root directory ----- #
current_dir = os.getcwd()
os.chdir(ic.dir_iraf)

from pyraf import iraf
from pyraf.iraf import gemini, gmos

os.chdir(current_dir)
iraf.chdir(current_dir)


# ----- Wavelength alignment ----- #
from astropy.io import fits
import pandas as pd

iraf.rv()
iraf.unlearn('rvidlines')
iraf.rv.observatory = ic.obs_site

# rvFlags = {'coordlist':'skylines.txt', 'ftype':'emission', 'nsum':1, 'threshold':7.,
#            'maxfeatures':13, 'fwidth':10., 'cradius':10., 'minsep':5.,
#            'logfile':'rvLog.txt', 'autowrite':'yes'}

for sci in iraf.type(ic.lst_sci, Stdout=1):
	sci = sci.strip()

	hd0 = fits.getheader('stxeqxbrg'+sci+'.fits', ext=0)
	hd2 = fits.getheader('stxeqxbrg'+sci+'.fits', ext=2)

	# fits.setval('stxeqxbrg'+sci+'.fits', 'DATE-OBS', value=hd0['DATE-OBS'], ext=5)
	# fits.setval('stxeqxbrg'+sci+'.fits', 'UT', value=hd0['UT'], ext=5)
	# fits.setval('stxeqxbrg'+sci+'.fits', 'RA', value=hd0['RA'], ext=5)
	# fits.setval('stxeqxbrg'+sci+'.fits', 'DEC', value=hd0['DEC'], ext=5)

	os.system('rm -rfv rvLog_'+sci+'.txt')

	# iraf.rvidlines('stxeqxbrg'+sci+'.fits[SKY]', **rvFlags)
	iraf.rvidlines('stxeqxbrg'+sci+'.fits[SKY]', coordlist='skylines.txt', ftype='emission',
                   nsum=1, threshold=7.0, maxfeatures=13, fwidth=10.0, cradius=10.0, minsep=5.0,
                   autowrite='yes', logfile='rvLog_'+sci+'.txt')

	f = open('rvLog_'+sci+'.txt','r')
	lines = f.readlines()
	bool_Zobs = pd.Series(lines).str.startswith('stxeqxbrg'+sci+'.fits[SKY]   1 : Zobs').values
	line_Zobs = np.array(lines)[bool_Zobs][0]
	Zobs = float(line_Zobs.split('Zobs     =')[1].split(',')[0])
	f.close()

	print('{0:.3f} --> {1:.3f}'.format(hd2['CRVAL1'], hd2['CRVAL1']+(-Zobs*hd0['CENTWAVE']*10.0)))
	nw_zero = hd2['CRVAL1']+(-Zobs*hd0['CENTWAVE']*10.0)

	iraf.unlearn('hedit')
	iraf.hedit.update = 'yes'
	iraf.hedit.verify = 'no'

	iraf.hedit('stxeqxbrg'+sci+'.fits[SCI]', 'CRVAL1', nw_zero)

