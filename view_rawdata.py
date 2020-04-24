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
dir_raw = 'Raw/'


# ----- Running DS9 for viewing raw data ----- #
os.chdir(dir_raw)
want_to_see = ['S20190228S0015',
               'S20190304S0134',
               'S20190304S0135']


for i in np.arange(len(want_to_see)):
	ds9_com = 'ds9 -multiframe '+want_to_see[i]+'.fits '
	ds9_opt = '-tile grid mode manual -tile grid layout 12 1 -scalemode zscale -scale lock yes -frame lock image &'
	os.system(ds9_com+ds9_opt)


os.chdir(current_dir)

