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
flat = np.loadtxt(ic.lst_flat, dtype=str)
if (flat.size > 1):
	raise ValueError("Please check if there is only one flat image for the standard star.")
flat0 = flat.item(0)

# Extract the arc
arc = np.loadtxt(ic.lst_arc, dtype=str)
if (arc.size > 1):
	raise ValueError("Please check if there is only one arc image for the standard star.")
arc0 = arc.item(0)

iraf.imdelete('g@'+ic.lst_arc)
iraf.imdelete('rg@'+ic.lst_arc)
iraf.imdelete('erg@'+ic.lst_arc)

iraf.gfreduce(arc0, rawpath=ic.rawdir, fl_extract='yes', recenter='no',
              trace='no', reference='erg'+flat0, fl_bias='no',
              fl_over='yes', fl_trim='yes', mdffile=ic.nmdf, mdfdir='./',
              slits=ic.cslit, fl_fluxcal='no', fl_gscrrej='no',
              fl_wavtran='no', fl_skysub='no', fl_inter='no')


# ----- Measure the wavelength solution ----- #
iraf.sleep(10.0)
iraf.gswavelength('erg'+arc0, fl_inter='yes',
                  nlost=10, ntarget=15, threshold=25,
                  coordlis='gmos$data/GCALcuar.dat')
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
