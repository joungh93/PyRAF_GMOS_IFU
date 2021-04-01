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
iraf.unlearn('gfdisplay')


# ---------- Trace reference ---------- #
flat = np.loadtxt(ic.lst_flat, dtype=str)
if (flat.size > 1):
	raise ValueError("Please check if there is only one flat image for the standard star.")
flat0 = flat.item(0)

iraf.imdelete('g@'+ic.lst_flat)
iraf.imdelete('rg@'+ic.lst_flat)
iraf.imdelete('erg@'+ic.lst_flat)

iraf.gfreduce(flat0, rawpath=ic.rawdir, fl_extract='yes', bias=ic.caldir+ic.procbias,
	          fl_over='yes', fl_trim='yes', mdffile=ic.nmdf, mdfdir='./',
	          slits=ic.cslit, line=pk_line, fl_fluxcal='no', fl_gscrrej='no',
	          fl_wavtran='no', fl_skysub='no', fl_inter='no', fl_vardq='yes')


# ----- Displaying the results ----- #
if (ic.nslit == 1):
	vkw = '1'
if (ic.nslit == 2):
	vkw = '*'
os.system('ds9 &')
iraf.sleep(5.0)
iraf.gfdisplay('erg'+flat0, 1, version=vkw)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
