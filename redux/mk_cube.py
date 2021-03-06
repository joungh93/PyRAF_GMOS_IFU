#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 16:29:33 2019

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

iraf.unlearn('gfcube')


# ----- Create the cubes ----- #
from astropy.io import fits

for sci in iraf.type(ic.lst_sci, Stdout=1):
	sci = sci.strip()
	iraf.imdelete('cstxeqxbrg'+sci+'_3D', verify='no')
	iraf.gfcube('cstxeqxbrg'+sci, outimage='cstxeqxbrg'+sci+'_3D',
		        fl_atmdisp='yes', fl_var='yes', fl_dq='yes')

	fits.open('cstxeqxbrg'+sci+'_3D.fits').info()

	os.system('ds9 cstxeqxbrg'+sci+'_3D.fits[sci] &')


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))

# # # Stop point #1 : pyfu is not working in the python script.
# # Import Error: cannot import name imagestats
# import sys
# sys.exit('Check the individual cubes!')


# # ----- Combine cubes (PyRAF scripts) ----- #
# cd ..
# iraf27
# pyraf
# pyfu
# cd redux3/

# import init_cfg as ic
# os.system('rm -rfv separatedcubes* '+ic.final_cube)
# pyfalign('*_3D*', llimit=200)
# pyfmosaic('*_3D*', 'separatedcubes', separate='yes')

# !ds9 &
# pyfmosaic('*_3D*', ic.final_cube, var='yes')
