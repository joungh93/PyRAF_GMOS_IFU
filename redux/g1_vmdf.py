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
dir_wav = sorted(glob.glob("w*"))
for i in np.arange(len(dir_wav)):
    if os.path.isdir(dir_wav[i]):
        pass
    else:
        raise ValueError("Please re-check the science directories.")

for i in np.arange(len(dir_wav)):
    dir_sci = sorted(glob.glob(dir_wav[i]+"/*"))

    for j in np.arange(1):#np.arange(len(dir_sci)):

        # Moving each science directory
        name_sci = dir_sci[j].split("/")[-1]
        print("Moving path for "+name_sci+"...")
        os.chdir(current_dir+"/"+dir_sci[j])
        iraf.chdir(current_dir+"/"+dir_sci[j])
        
        os.system("rm -rfv "+ic.dir_db+" *.fits tmp*")
        os.system("cp -rpv "+ic.dir_std+ic.nmdf+" .")

        flat = np.loadtxt(ic.lst_flat, dtype=str)
        flat0 = flat.item(0)

        iraf.imdelete('g@'+ic.lst_flat)
        iraf.imdelete('rg@'+ic.lst_flat)
        iraf.gfreduce(flat0, rawpath=ic.rawdir, fl_extract='no', bias=ic.caldir+ic.procbias,
                      fl_over='yes', fl_trim='yes', mdffile=ic.nmdf, mdfdir='./',
                      slits=ic.cslit, line=ic.pk_line, fl_fluxcal='no', fl_gscrrej='no',
                      fl_wavtran='no', fl_skysub='no', fl_inter='no', fl_vardq='no')

        stdflat = np.loadtxt(ic.dir_std+ic.lst_stdflat, dtype=str)
        stdflat = stdflat.item(0)

        # if (glob.glob(ic.dir_db) == []):
        #     os.system("mkdir "+ic.dir_db)

        # com_cpdb = "cp -rpv "
        # com_cpdb += ic.dir_std+ic.dir_db+"aperg"+stdflat+"_1 "
        # com_cpdb += ic.dir_db+"aperg"+flat0+"_1"
        # os.system(com_cpdb)
        # if (ic.nslit == 2):
        #     com_cpdb2 = "cp -rpv "
        #     com_cpdb2 += ic.dir_std+ic.dir_db+"aperg"+stdflat+"_2 "
        #     com_cpdb2 += ic.dir_db+"aperg"+flat0+"_2"
        #     os.system(com_cpdb2) 

        if (i == 1):           
            iraf.imdelete('erg@'+ic.lst_flat)
            iraf.gfextract('rg'+flat0, fl_inter='yes', line=ic.pk_line, exslits=ic.eslit)

        # Coming back to current path
        os.chdir(current_dir)
        iraf.chdir(current_dir)     



'''

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


# Writing aperture file
if (ic.nslit == 1):
    apfile = ['aperg'+flat0+'_1']
if (ic.nslit == 2):
    apfile = ['aperg'+flat0+'_1', 'aperg'+flat0+'_2']
for i in np.arange(len(apfile)):
    os.system('cp -rpv '+ic.dir_db+apfile[i]+' '+ic.dir_db+apfile[i]+'_old')

'''

# Printing the running time
print('--- %.4f seconds ---' %(time.time()-start_time))



