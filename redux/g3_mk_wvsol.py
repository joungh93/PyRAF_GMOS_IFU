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

iraf.unlearn('gfreduce')
iraf.unlearn('gswavelength')


# ---------- Wavelength solution ---------- #
for d in ic.dir_wav:
    dir_sci = sorted(glob.glob(d+"/*"))

    for j in np.arange(len(dir_sci)):

        # Moving each science directory
        name_sci = dir_sci[j].split("/")[-1]
        print("Moving path for "+name_sci+"...")
        os.chdir(current_dir+"/"+dir_sci[j])
        iraf.chdir(current_dir+"/"+dir_sci[j])

        # FLAT
        flat = np.loadtxt(ic.lst_flat, dtype=str)
        flat0 = flat.item(0)

        # ARC
        arc = np.loadtxt(ic.lst_arc, dtype=str)
        arc0 = arc.item(0)

        # Tracing arc
        iraf.imdelete('g@'+ic.lst_arc)
        iraf.imdelete('rg@'+ic.lst_arc)
        iraf.imdelete('erg@'+ic.lst_arc)
        iraf.gfreduce(arc0, rawpath=ic.rawdir, fl_extract='yes', recenter='no',
                      trace='no', reference='erg'+flat0, fl_bias='no',
                      fl_over='yes', fl_trim='yes', mdffile=ic.nmdf, mdfdir='./',
                      slits=ic.cslit, fl_fluxcal='no', fl_gscrrej='no',
                      fl_wavtran='no', fl_skysub='no', fl_inter='no')

        # Measure the wavelength solution
        iraf.sleep(2.0)
        iraf.gswavelength('erg'+arc0, fl_inter='yes',
                          # nlost=10, ntarget=15, threshold=25,
                          nlost=20, ntarget=30, threshold=0,
                          coordlis='gmos$data/GCALcuar.dat')        

        # Coming back to current path
        os.chdir(current_dir)
        iraf.chdir(current_dir)  
        
'''
----- Interactive task after gswavelength -----
Examine identifications interactively? (Enter)
(IRAF graphics of spectrum displaying...)

"The spectrum window"
    - "w" + "e" (left bottom) + "e" (right top) : zoom-in
    - "w" + "a" : zoom-out
    - "d" :  delete the line
    - "m" : mark the line
    - "f" : jump to the parabola window
    - "q" : quitting the interactive task

"The parabola window"
    - "d" : jump to the spectrum window
    - "f" : fit the line again
    - "q" : return to the spectrum window

For the two-slit mode, you have to do the manual check twice.

Fit dispersion function interactively? (no|yes|NO|YES) ('NO'): Enter
Output file : erg[ARC].fits, database/aperg[ARC]_[1,2], database/iderg[ARC]_[1,2]
'''


# Printing the running time
print('--- %.4f seconds ---' %(time.time()-start_time))
