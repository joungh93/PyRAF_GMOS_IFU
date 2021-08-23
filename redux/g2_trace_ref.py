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
for d in ic.dir_wav:
    dir_sci = sorted(glob.glob(d+"/*"))

    for j in np.arange(len(dir_sci)):

        # Moving each science directory
        name_sci = dir_sci[j].split("/")[-1]
        print("Moving path for "+name_sci+"...")
        os.chdir(current_dir+"/"+dir_sci[j])
        iraf.chdir(current_dir+"/"+dir_sci[j])

        # Trace reference from flat data
        flat = np.loadtxt(ic.lst_flat, dtype=str)
        flat0 = flat.item(0)

        iraf.imdelete('g@'+ic.lst_flat)
        iraf.imdelete('rg@'+ic.lst_flat)
        iraf.imdelete('erg@'+ic.lst_flat)
        iraf.gfreduce(flat0, rawpath=ic.rawdir, fl_extract='yes', bias=ic.caldir+ic.procbias,
                      fl_over='yes', fl_trim='yes', mdffile=ic.nmdf, mdfdir='./',
                      slits=ic.cslit, line=ic.pk_line, fl_fluxcal='no', fl_gscrrej='no',
                      fl_wavtran='no', fl_skysub='no', fl_inter='no', fl_vardq='yes')

        # Coming back to current path
        os.chdir(current_dir)
        iraf.chdir(current_dir)     


# Printing the running time
print('--- %.4f seconds ---' %(time.time()-start_time))
