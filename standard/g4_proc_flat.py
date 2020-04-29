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

# Model and remove the light

iraf.imdelete('brg@'+ic.lst_flat)

# os.system('ds9 &')
# iraf.sleep(5.0)
for flat in iraf.type(ic.lst_flat, Stdout=1):
    flat = flat.strip()
    blkmsk = 'newblkmask_S20190228S0015_hdr01'
    iraf.gfscatsub('rg'+flat, blkmsk, outimage='', prefix='b',
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
iraf.imdelete('eqbrg@'+ic.lst_flat)

flatref = iraf.type(ic.lst_flat, Stdout=1)[0]
arc0 = iraf.type(ic.lst_arc, Stdout=1)[0]

iraf.imdelete('qbrg@'+ic.lst_flat)
iraf.imdelete('eqbrg@'+ic.lst_flat)
iraf.gfreduce('brg@'+ic.lst_flat, recenter='no', reference='erg'+flatref,
              fl_extract='yes', fl_qecorr='yes', qe_refim='erg'+arc0,
              fl_addmdf='no', fl_bias='no', fl_over='no', fl_trim='no',
              mdffile=ic.nmdf, mdfdir='./',
              slits='both', line=1400, fl_fluxcal='no', fl_gscrrej='no',
              fl_wavtran='no', fl_skysub='no', fl_inter='no',
              fl_vardq='yes')

os.system('ds9 &')
iraf.sleep(5.0)
for flat in iraf.type(ic.lst_flat, Stdout=1):
	flat = flat.strip()
	iraf.gfdisplay('eqbrg'+flat, 1)


# ---------- Response function ---------- #
for flat in iraf.type(ic.lst_flat, Stdout=1):
    flat = flat.strip()
    iraf.imdelete(flat+'_resp')
    iraf.gfresponse('eqbrg'+flat, outimage=flat+'_resp', sky='', 
    	              order=45, func='spline3', sample='*', 
    	              fl_fit='yes', fl_inter='no')

# os.system('ds9 &')
# iraf.sleep(5.0)
for flat in iraf.type(ic.lst_flat, Stdout=1):
    flat = flat.strip()
    iraf.gfdisplay(flat+'_resp', 1)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
