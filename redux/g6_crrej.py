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

iraf.unlearn('gemcrspec')


# ---------- Cosmic ray rejection ---------- #
for d in ic.dir_wav:
    dir_sci = sorted(glob.glob(d+"/*"))

    for j in np.arange(len(dir_sci)):

        # Moving each science directory
        name_sci = dir_sci[j].split("/")[-1]
        print("Moving path for "+name_sci+"...")
        os.chdir(current_dir+"/"+dir_sci[j])
        iraf.chdir(current_dir+"/"+dir_sci[j])

        # Running gemcrspec
        sci = np.loadtxt(ic.lst_sci, dtype=str)
        sci0 = sci.item(0)

        iraf.imdelete('xbrg@'+ic.lst_sci)
        iraf.gemcrspec('brg'+sci0, 'xbrg'+sci0, logfile='crrej.log',
                       key_gain='GAIN', key_ron='RDNOISE', xorder=9,
                       yorder=-1, sigclip=4.5, sigfrac=0.5, objlim=1.,
                       niter=4, verbose='yes', fl_vardq='yes')

        # Coming back to current path
        os.chdir(current_dir)
        iraf.chdir(current_dir)          


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
