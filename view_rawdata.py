#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 2 14:46:15 2020

@author: jlee
"""


import numpy as np
import glob, os
from astropy.io import fits


current_dir = os.getcwd()
dir_raw = 'raw/'
use_IFU2 = False


# ----- Running DS9 for viewing raw data ----- #
os.chdir(dir_raw)
want_to_see = ['N20141015S0456',
               'N20141217S0199',
               'N20141217S0200',
               'N20141217S0201',
               'N20141221S0207',
               'N20141222S0521']

if use_IFU2:
	n_grid = 12
else:
	n_grid = 6

for i in np.arange(len(want_to_see)):
	ds9_com = 'ds9 -multiframe '+want_to_see[i]+'.fits '
	ds9_opt = '-tile grid mode manual -tile grid layout '+str(n_grid)+' 1 '
	ds9_opt += '-scalemode zscale -scale lock yes -frame lock image &'
	os.system(ds9_com+ds9_opt)


os.chdir(current_dir)

