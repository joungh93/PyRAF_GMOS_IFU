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

iraf.unlearn('gscalibrate')


# ----- Spectrophotometric calibration ----- #
for d in ic.dir_wav:
    dir_sci = sorted(glob.glob(d+"/*"))

    for j in np.arange(len(dir_sci)):

        # Moving each science directory
        name_sci = dir_sci[j].split("/")[-1]
        print("Moving path for "+name_sci+"...")
        os.chdir(current_dir+"/"+dir_sci[j])
        iraf.chdir(current_dir+"/"+dir_sci[j])

        # SCI
        sci = np.loadtxt(ic.lst_sci, dtype=str)
        sci0 = sci.item(0)

        # Flux calibration
        iraf.imdelete('cstxeqxbrg@'+ic.lst_sci, verify='no')
        iraf.gscalibrate('stxeqxbrg'+sci0, sfunction=ic.caldir+ic.sensfunc,
                         obs=ic.obs_site, extinction=ic.extinction,
                         fl_ext='yes', fl_vardq='yes')        

        # Coming back to current path
        os.chdir(current_dir)
        iraf.chdir(current_dir)  


# os.system('ds9 &')
# iraf.sleep(5.0)
# for sci in iraf.type(ic.lst_sci, Stdout=1):
#   iraf.gfdisplay('cstxeqxbrg'+sci, 1, version=1)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
