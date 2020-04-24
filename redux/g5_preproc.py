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
iraf.unlearn('gfscatsub')


# ---------- Pre-processing of the science frames ---------- #

# MDF, bias, and overscan
iraf.imdelete('g@'+ic.lst_sci)
iraf.imdelete('rg@'+ic.lst_sci)

iraf.gfreduce('@'+ic.lst_sci, rawpath=ic.rawdir, fl_extract='no',
	            bias=ic.bias, fl_over='yes', fl_trim='yes', mdffile=ic.nmdf, mdfdir='./',
              slits='both', line=1400, fl_fluxcal='no', fl_gscrrej='no',
              fl_wavtran='no', fl_skysub='no', fl_vardq='yes', fl_inter='no')


# Scattered light
blkmsk_ref = 'newblkmask_S20190228S0015_hdr01'

iraf.imdelete('brg@'+ic.lst_sci)

for sci in iraf.type(ic.lst_sci, Stdout=1):
    sci = sci.strip()
    iraf.gfscatsub('rg'+sci, blkmsk_ref, outimage='', prefix='b',
                   xorder='3,3,3,3,3,3,3,3,3,3,3,3',
                   yorder='3,3,3,3,3,3,3,3,3,3,3,3',
                   cross='yes', fl_inter='no')

# os.system('ds9 &')
# iraf.sleep(5.0)
# for sci in iraf.type(ic.lst_sci, Stdout=1):
#     sci = sci.strip()
#     for i in range(12):
#         iraf.imexamine('brg'+sci+'[sci,'+str(i+1)+']', 1)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
