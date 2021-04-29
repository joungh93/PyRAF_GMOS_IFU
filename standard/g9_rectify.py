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
from astropy.io import fits


# ----- Importing IRAF from the root directory ----- #
current_dir = os.getcwd()
os.chdir(ic.dir_iraf)

from pyraf import iraf
from pyraf.iraf import gemini, gmos

os.chdir(current_dir)
iraf.chdir(current_dir)

iraf.unlearn('gftransform')
iraf.unlearn('gfskysub')


# ---------- Rectify the spectra ---------- #

# Angstroms per pixel
arc = np.loadtxt(ic.lst_arc, dtype=str)
if (arc.size > 1):
    raise ValueError("Please check if there is only one arc image for the standard star.")
arc0 = arc.item(0)

std = np.loadtxt(ic.lst_std, dtype=str)
if (std.size > 1):
    raise ValueError("Please check if there is only one image for the standard star.")
std0 = std.item(0)

iraf.imdelete('txeqxbrg@'+ic.lst_std, verify='no')
iraf.gftransform('xeqxbrg'+std0, wavtraname='erg'+arc0, fl_vardq='yes')


# ---------- Sky subtraction ---------- #
iraf.imdelete('stxeqxbrg@'+ic.lst_std, verify='no')
iraf.gfskysub('txeqxbrg'+std0, fl_inter='no', combine='median', sepslits='yes')

ds9_comm = "ds9 -scalemode zscale -scale lock yes -frame lock image "
os.system(ds9_comm+"txeqxbrg"+std0+".fits[2] stxeqxbrg"+std0+".fits[2] &")


# Printing the running time
print('--- %.4f seconds ---' %(time.time()-start_time))
