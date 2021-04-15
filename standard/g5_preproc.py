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

iraf.unlearn('gfreduce')
iraf.unlearn('gfscatsub')


# ---------- Pre-processing of the science frames ---------- #

# MDF, bias, and overscan
iraf.imdelete('g@'+ic.lst_std)
iraf.imdelete('rg@'+ic.lst_std)

iraf.gfreduce('@'+ic.lst_std, rawpath=ic.rawdir, fl_extract='no',
	          bias=ic.caldir+ic.procbias, fl_over='yes', fl_trim='yes',
	          mdffile=ic.nmdf, mdfdir='./',
              slits=ic.cslit, line=pk_line, fl_fluxcal='no', fl_gscrrej='no',
              fl_wavtran='no', fl_skysub='no', fl_vardq='yes', fl_inter='no')


# Scattered light
blkmsk = np.loadtxt("blkmask_name.txt", dtype=str).item(0)
blkmsk0 = blkmsk

iraf.imdelete('brg@'+ic.lst_std)

std = np.loadtxt(ic.lst_std, dtype=str)
if (std.size > 1):
    raise ValueError("Please check if there is only one image for the standard star.")
std0 = std.item(0)

iraf.gfscatsub('rg'+std0, blkmsk0, outimage='', prefix='b',
               xorder='3,3,3,3,3,3,3,3,3,3,3,3',
               yorder='3,3,3,3,3,3,3,3,3,3,3,3',
               cross='yes', fl_inter='no')

# os.system('ds9 &')
# iraf.sleep(5.0)
# for std in iraf.type(ic.lst_std, Stdout=1):
#     std = std.strip()
#     for i in range(12):
#         iraf.imexamine('brg'+std+'[sci,'+str(i+1)+']', 1)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
