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

iraf.unlearn('gftransform')
iraf.unlearn('gfskysub')


# ---------- Rectify the spectra ---------- #
for d in ic.dir_wav:
    dir_sci = sorted(glob.glob(d+"/*"))

    for j in np.arange(len(dir_sci)):

        # Moving each science directory
        name_sci = dir_sci[j].split("/")[-1]
        print("Moving path for "+name_sci+"...")
        os.chdir(current_dir+"/"+dir_sci[j])
        iraf.chdir(current_dir+"/"+dir_sci[j])

        # SCI & ARC
        sci = np.loadtxt(ic.lst_sci, dtype=str)
        sci0 = sci.item(0)

        arc = np.loadtxt(ic.lst_arc, dtype=str)
        arc0 = arc.item(0)

        # Rectify the spectra
        iraf.imdelete('txeqxbrg@'+ic.lst_sci, verify='no')
        iraf.gftransform('xeqxbrg'+sci0, wavtraname='erg'+arc0, fl_vardq='yes')

        # Sky subtraction
        iraf.imdelete('stxeqxbrg@'+ic.lst_sci, verify='no')
        iraf.gfskysub('txeqxbrg'+sci0, fl_inter='no', combine='median', sepslits='yes')

        ds9_comm = "ds9 -scalemode zscale -scale lock yes -frame lock image "
        os.system(ds9_comm+"txeqxbrg"+sci0+".fits[2] stxeqxbrg"+sci0+".fits[2] &")

        # Coming back to current path
        os.chdir(current_dir)
        iraf.chdir(current_dir)  


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
