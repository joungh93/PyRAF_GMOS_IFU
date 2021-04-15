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

iraf.unlearn('gfscatsub')
iraf.unlearn('gfreduce')
iraf.unlearn('gfdisplay')
iraf.unlearn('gfresponse')


# ---------- Reduce the lamp flat ---------- #
flat = np.loadtxt(ic.lst_flat, dtype=str)
if (flat.size > 1):
    raise ValueError("Please check if there is only one flat image for the standard star.")
flat0 = flat.item(0)

arc = np.loadtxt(ic.lst_arc, dtype=str)
if (arc.size > 1):
    raise ValueError("Please check if there is only one arc image for the standard star.")
arc0 = arc.item(0)

# Model and remove the light
iraf.imdelete('brg@'+ic.lst_flat)

# os.system('ds9 &')
# iraf.sleep(5.0)
blkmsk = np.loadtxt("blkmask_name.txt", dtype=str).item(0)
blkmsk0 = blkmsk
iraf.gfscatsub('rg'+flat0, blkmsk0, outimage='', prefix='b',
               xorder='3,3,3,3,3,3,3,3,3,3,3,3',
               yorder='3,3,3,3,3,3,3,3,3,3,3,3', 
               cross='yes', fl_inter='no')

# os.system('ds9 &')
# iraf.sleep(5.0)
# for flat in iraf.type(ic.lst_flat, Stdout=1):
#     flat = flat.strip()
#     for i in np.arange(12):
#         iraf.imexamine('brg'+flat+'[sci,'+str(i+1)+']', 1)


# QE correction and extract
iraf.imdelete('qbrg@'+ic.lst_flat)
iraf.imdelete('eqbrg@'+ic.lst_flat)
iraf.gfreduce('brg@'+ic.lst_flat, recenter='no', reference='erg'+flat0,
              fl_extract='yes', fl_qecorr='yes', qe_refim='erg'+arc0,
              fl_addmdf='no', fl_bias='no', fl_over='no', fl_trim='no',
              mdffile=ic.nmdf, mdfdir='./',
              slits=ic.cslit, line=pk_line, fl_fluxcal='no', fl_gscrrej='no',
              fl_wavtran='no', fl_skysub='no', fl_inter='no', fl_vardq='yes')

if (ic.nslit == 1):
    vkw = '1'
if (ic.nslit == 2):
    vkw = '*'
os.system('ds9 &')
iraf.sleep(5.0)
iraf.gfdisplay('eqbrg'+flat0, 1, version=vkw)


# ---------- Response function ---------- #
for flat in iraf.type(ic.lst_flat, Stdout=1):
    flat = flat.strip()
    iraf.imdelete(flat0+'_resp')
    iraf.gfresponse('eqbrg'+flat0, outimage=flat0+'_resp', sky='', 
    	            order=45, func='spline3', sample='*', 
    	            fl_fit='yes', fl_inter='no')

# os.system('ds9 &')
# iraf.sleep(5.0)
iraf.gfdisplay(flat0+'_resp', 1, version=vkw)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
