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

# Extract the arc
flatref = iraf.type(ic.lst_flat, Stdout=1)[0].strip()

iraf.imdelete('g@'+ic.lst_arc)
iraf.imdelete('rg@'+ic.lst_arc)
iraf.imdelete('erg@'+ic.lst_arc)

for arc in iraf.type(ic.lst_arc, Stdout=1):
	arc = arc.strip()
	iraf.gfreduce(arc, rawpath=ic.rawdir, fl_extract='yes', recenter='no',
	              trace='no', reference='erg'+flatref, fl_bias='no',
	              fl_over='yes', fl_trim='yes', mdffile=ic.nmdf, mdfdir='./',
	              slits='both', fl_fluxcal='no', fl_gscrrej='no',
	              fl_wavtran='no', fl_skysub='no', fl_inter='no')

iraf.sleep(10.0)

# ----- Measure the wavelength solution ----- #
for arc in iraf.type(ic.lst_arc, Stdout=1):
	arc = arc.strip()
	iraf.gswavelength('erg'+arc, fl_inter='yes',
	                  nlost=20, ntarget=15, threshold=25,
	                  coordlis='gmos$data/GCALcuar.dat')


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
