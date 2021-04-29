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

iraf.imdelete('txeqxbrg'+std0, verify='no')
iraf.gftransform('xeqxbrg'+std0, wavtraname='erg'+arc0, fl_vardq='yes')


# ---------- Sky subtraction ---------- #
iraf.imdelete('stxeqxbrg@'+ic.lst_std, verify='no')
iraf.gfskysub('txeqxbrg'+std0, fl_inter='no',
              combine='median', sepslits='yes')

os.system('ds9 &')
iraf.sleep(5.0)
for std in iraf.type(ic.lst_std, Stdout=1):
    std = std.strip()
    iraf.display('txeqxbrg'+std+'[sci,1]', 1)
    iraf.display('stxeqxbrg'+std+'[sci,1]', 2)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
