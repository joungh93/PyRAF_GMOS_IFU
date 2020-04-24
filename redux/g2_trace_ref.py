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
iraf.unlearn('gfdisplay')


# ---------- Trace reference ---------- #
iraf.imdelete('g@'+ic.lst_flat)
iraf.imdelete('rg@'+ic.lst_flat)
iraf.imdelete('erg@'+ic.lst_flat)

for flat in iraf.type(ic.lst_flat, Stdout=1):
	flat = flat.strip()
	iraf.gfreduce(flat, rawpath=ic.rawdir, fl_extract='yes', bias=ic.bias,
                  fl_over='yes', fl_trim='yes', mdffile=ic.nmdf, mdfdir='./',
                  slits='both', line=1400, fl_fluxcal='no', fl_gscrrej='no',
                  fl_wavtran='no', fl_skysub='no', fl_inter='no', fl_vardq='yes')

os.system('ds9 &')
iraf.sleep(5.0)
for flat in iraf.type(ic.lst_flat, Stdout=1):
    flat = flat.strip()
    iraf.gfdisplay('erg'+flat, 1)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
