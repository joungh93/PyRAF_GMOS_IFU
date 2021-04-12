#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 17:21:30 2019

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import glob, os
import g0_init_cfg as ic


# ----- Line number (to be revised!) ----- #
pk_line = 1400
'''
Line for finding peaks (gfreduce)
Line/column for finding apertures (gfextract)
'''
# ---------------------------------------- #


# ----- Importing IRAF from the root directory ----- #
current_dir = os.getcwd()
os.chdir(ic.dir_iraf)

from pyraf import iraf
from pyraf.iraf import gemini, gmos

os.chdir(current_dir)
iraf.chdir(current_dir)

iraf.unlearn('gfreduce')
iraf.unlearn('gfextract')


# ---------- Verifying the MDF ---------- #
os.system('rm -rfv '+ic.dir_db+' *.fits tmp*')

# Copy the MDF
iraf.copy('gmos$data/'+ic.mdf, '.', verbose='no')

# Extract a flat
flat = np.loadtxt(ic.lst_flat, dtype=str)
if (flat.size > 1):
	raise ValueError("Please check if there is only one flat image for the standard star.")
flat0 = flat.item(0)

iraf.imdelete('g@'+ic.lst_flat)
iraf.imdelete('rg@'+ic.lst_flat)
iraf.gfreduce(flat0, rawpath=ic.rawdir, fl_extract='no', bias=ic.caldir+ic.procbias,
	          fl_over='yes', fl_trim='yes', mdffile=ic.mdf, mdfdir='./',
	          slits=ic.cslit, line=pk_line, fl_fluxcal='no', fl_gscrrej='no',
	          fl_wavtran='no', fl_skysub='no', fl_inter='no', fl_vardq='no')

iraf.imdelete('erg@'+ic.lst_flat)
iraf.gfextract('rg'+flat0, fl_inter='yes', line=pk_line, exslits=ic.eslit)
'''
----- Interactive task after gfextract -----
Extracting slit 1
Find apertures for erg[FLAT]_1? ('yes')
Edit apertures for erg[FLAT]_1? ('yes')
(IRAF graphics displaying... please check the fibers visually.)
- "w" + "e" (left bottom) + "e" (right top) : zoom-in
- "w" + "a" : zoom-out
- "q" : quitting the interactive task

Trace apertures for erg[FLAT]_1? ('yes')
Fit traced positions for erg[FLAT]_1 interactively? ('NO')
Write apertures for erg[FLAT]_1 to database ('yes')
Extract aperture spectra for erg[FLAT]_1? ('yes' (IFU-2 slit) / 'NO' (IFU-1 slit))
--> For the IFU-1 slit, this is the end.
Review extracted spectra from erg[FLAT]_1? ('NO')

Extracting slit 2
Find apertures for erg[FLAT]_2? ('yes')
(...repeating the tasks for slit 2...)
Extract aperture spectra for erg[FLAT]_1? ('NO')

--> 'GFEXTRACT exit status: error' message will appear, but it is not a fault.
'''

# Writing aperture file
if (ic.nslit == 1):
	apfile = ['aperg'+flat0+'_1']
if (ic.nslit == 2):
	apfile = ['aperg'+flat0+'_1', 'aperg'+flat0+'_2']
for i in np.arange(len(apfile)):
	os.system('cp -rpv '+ic.dir_db+apfile[i]+' '+ic.dir_db+apfile[i]+'_old')


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
