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

iraf.unlearn('gfscatsub')
iraf.unlearn('gfreduce')
iraf.unlearn('gfdisplay')
iraf.unlearn('gfresponse')


# ---------- Reduce the lamp flat ---------- #
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

        # Model and remove the light
        iraf.imdelete('brg@'+ic.lst_flat)
        blkmsk = np.loadtxt("blkmask_name.txt", dtype=str)
        blkmsk0 = blkmsk.item(0)
        iraf.gfscatsub('rg'+flat0, blkmsk0, outimage='', prefix='b',
                       xorder='3,3,3,3,3,3,3,3,3,3,3,3',
                       yorder='3,3,3,3,3,3,3,3,3,3,3,3', 
                       cross='yes', fl_inter='no')

        # QE correction and extract
        iraf.imdelete('qbrg@'+ic.lst_flat)
        iraf.imdelete('eqbrg@'+ic.lst_flat)
        iraf.gfreduce('brg@'+ic.lst_flat, recenter='no', reference='erg'+flat0,
                      fl_extract='yes', fl_qecorr='yes', qe_refim='erg'+arc0,
                      fl_addmdf='no', fl_bias='no', fl_over='no', fl_trim='no',
                      mdffile=ic.nmdf, mdfdir='./',
                      slits=ic.cslit, line=ic.pk_line, fl_fluxcal='no', fl_gscrrej='no',
                      fl_wavtran='no', fl_skysub='no', fl_inter='no', fl_vardq='yes')

        # Response function
        iraf.imdelete(flat0+'_resp')
        iraf.gfresponse('eqbrg'+flat0, outimage=flat0+'_resp', sky='', 
                        order=45, func='spline3', sample='*', 
                        fl_fit='yes', fl_inter='no')

        # Coming back to current path
        os.chdir(current_dir)
        iraf.chdir(current_dir)  


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
