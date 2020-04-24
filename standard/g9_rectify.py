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

iraf.unlearn('gftransform')
iraf.unlearn('gfskysub')


# ---------- Rectify the spectra ---------- #
from astropy.io import fits

# Angstroms per pixel
arc0 = iraf.type(ic.lst_arc, Stdout=1)[0].strip()

for std in iraf.type(ic.lst_std, Stdout=1):
    std = std.strip()
    iraf.imdelete('txeqxbrg'+std, verify='no')
    iraf.gftransform('xeqxbrg'+std, wavtraname='erg'+arc0, fl_vardq='no')
    fits.open('txeqxbrg'+std+'.fits').info()
    dat, hdr = fits.getdata('txeqxbrg'+std+'.fits', ext=2, header=True)
    dw = float(hdr['CD1_1'])
    print('dw : {0:f}'.format(dw))


# # Stop point #1 
# import sys
# sys.exit("Please check 'dw'.")


# Rectify
dw = 1.93

iraf.imdelete('txeqxbrg@'+ic.lst_std)

for std in iraf.type(ic.lst_std, Stdout=1):
    std = std.strip()
    iraf.gftransform('xeqxbrg'+std, wavtraname='erg'+arc0, dw=dw, fl_vardq='yes')

# os.system('ds9 &')
# iraf.sleep(5.0)
# for std in iraf.type(ic.lst_std, Stdout=1):
#     std = std.strip()
#     iraf.display('txeqxbrg'+std+'.fits[sci,1]', 1)


# ---------- Sky subtraction ---------- #
iraf.imdelete('stxeqxbrg@'+ic.lst_std, verify='no')

for std in iraf.type(ic.lst_std, Stdout=1):
    std = std.strip()
    iraf.gfskysub('txeqxbrg'+std, fl_inter='no')

os.system('ds9 &')
iraf.sleep(5.0)
for std in iraf.type(ic.lst_std, Stdout=1):
    std = std.strip()
    iraf.display('txeqxbrg'+std+'[sci,1]', 1)
    iraf.display('stxeqxbrg'+std+'[sci,1]', 2)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))