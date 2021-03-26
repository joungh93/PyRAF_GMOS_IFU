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
os.system('rm -rfv '+ic.dir_db+' *.fits tmp*')

# Copy the MDF
iraf.copy('gmos$data/'+ic.mdf, '.', verbose='no')

# Extract a flat
iraf.imdelete('g@'+ic.lst_flat)
iraf.imdelete('rg@'+ic.lst_flat)

flat0 = iraf.type(ic.lst_flat, Stdout=1)[0]
iraf.gfreduce(flat0, rawpath=ic.rawdir, fl_extract='no', bias=ic.caldir+ic.procbias,
	          fl_over='yes', fl_trim='yes', mdffile=ic.mdf, mdfdir='./',
	          slits='both', line=1400, fl_fluxcal='no', fl_gscrrej='no',
	          fl_wavtran='no', fl_skysub='no', fl_inter='no', fl_vardq='no')

iraf.imdelete('erg@'+ic.lst_flat)
iraf.gfextract('rg'+flat0, fl_inter='yes', line=1400)


# Writing aperture file
apfile = ['aperg'+flat0+'_1', 'aperg'+flat0+'_2']
for i in np.arange(len(apfile)):
	os.system('cp -rpv '+ic.dir_db+apfile[i]+' '+ic.dir_db+apfile[i]+'_old')


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
