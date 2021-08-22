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
for d in ic.dir_wav:
    dir_sci = sorted(glob.glob(d+"/*"))

    for j in np.arange(len(dir_sci)):

        # Moving each science directory
        name_sci = dir_sci[j].split("/")[-1]
        print("Moving path for "+name_sci+"...")
        os.chdir(current_dir+"/"+dir_sci[j])
        iraf.chdir(current_dir+"/"+dir_sci[j])
        
        # Copying the MDF from standard star
        os.system("rm -rfv "+ic.dir_db+" *.fits tmp*")
        os.system("cp -rpv "+ic.dir_std+ic.nmdf+" .")

        # Verify the MDF from flat data
        flat = np.loadtxt(ic.lst_flat, dtype=str)
        flat0 = flat.item(0)

        iraf.imdelete('g@'+ic.lst_flat)
        iraf.imdelete('rg@'+ic.lst_flat)
        iraf.gfreduce(flat0, rawpath=ic.rawdir, fl_extract='no', bias=ic.caldir+ic.procbias,
                      fl_over='yes', fl_trim='yes', mdffile=ic.nmdf, mdfdir='./',
                      slits=ic.cslit, line=ic.pk_line, fl_fluxcal='no', fl_gscrrej='no',
                      fl_wavtran='no', fl_skysub='no', fl_inter='no', fl_vardq='no')

        # Interative tasks for the first science data for each central wavelength
        if (j == 0):
            fl_inter = 'yes'
        else:
            fl_inter = 'no' 
        iraf.imdelete('erg@'+ic.lst_flat)
        iraf.gfextract('rg'+flat0, fl_inter=fl_inter, line=ic.pk_line, exslits=ic.eslit)

        # Writing aperture file
        if (ic.nslit == 1):
            apfile = ['aperg'+flat0+'_1']
        if (ic.nslit == 2):
            apfile = ['aperg'+flat0+'_1', 'aperg'+flat0+'_2']
        for k in np.arange(len(apfile)):
            os.system('cp -rpv '+ic.dir_db+apfile[k]+' '+ic.dir_db+apfile[k]+'_old')

        # Coming back to current path
        os.chdir(current_dir)
        iraf.chdir(current_dir)     
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


# Printing the running time
print('--- %.4f seconds ---' %(time.time()-start_time))
