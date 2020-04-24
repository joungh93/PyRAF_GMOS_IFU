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

iraf.unlearn('gqecorr')
iraf.unlearn('gfextract')
iraf.unlearn('gfdisplay')


# ---------- QE correction and extraction of the science ---------- #
flat0 = iraf.type(ic.lst_flat, Stdout=1)[0].strip()
response = flat0+'_resp'
ref_flat0 = 'eqbrg'+flat0

arc0 = iraf.type(ic.lst_arc, Stdout=1)[0].strip()

iraf.imdelete('qxbrg@'+ic.lst_std)
iraf.imdelete('eqxbrg@'+ic.lst_std)

for std in iraf.type(ic.lst_std, Stdout=1):
    std = std.strip()
    iraf.gqecorr('xbrg'+std, refimage='erg'+arc0, fl_correct='yes',
                 fl_vardq='yes', verbose='yes')
    iraf.gfextract('qxbrg'+std, response=response, recenter='no',
                   trace='no', reference=ref_flat0, weights='none',
                   fl_vardq='yes', line=1400)

os.system('ds9 &')
iraf.sleep(5.0)
for std in iraf.type(ic.lst_std, Stdout=1):
	std = std.strip()
	iraf.gfdisplay('eqxbrg'+std, 1)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))