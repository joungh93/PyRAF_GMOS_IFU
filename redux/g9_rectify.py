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

for sci in iraf.type(ic.lst_sci, Stdout=1):
    sci = sci.strip()
    iraf.imdelete('txeqxbrg'+sci, verify='no')
    iraf.gftransform('xeqxbrg'+sci, wavtraname='erg'+arc0, fl_vardq='no')
    fits.open('txeqxbrg'+sci+'.fits').info()
    dat, hdr = fits.getdata('txeqxbrg'+sci+'.fits', ext=2, header=True)
    dw = float(hdr['CD1_1'])
    print('dw : {0:f}'.format(dw))


# # Stop point #1 
# import sys
# sys.exit("Please check 'dw'.")


# Rectify
dw = 1.93

iraf.imdelete('txeqxbrg@'+ic.lst_sci)

for sci in iraf.type(ic.lst_sci, Stdout=1):
    sci = sci.strip()
    iraf.gftransform('xeqxbrg'+sci, wavtraname='erg'+arc0, dw=dw, fl_vardq='yes')

# os.system('ds9 &')
# iraf.sleep(5.0)
# for sci in iraf.type(ic.lst_sci, Stdout=1):
#     sci = sci.strip()
#     iraf.display('txeqxbrg'+sci+'.fits[sci,1]', 1)


# ---------- Sky subtraction ---------- #
iraf.imdelete('stxeqxbrg@'+ic.lst_sci, verify='no')

for sci in iraf.type(ic.lst_sci, Stdout=1):
    sci = sci.strip()
    iraf.gfskysub('txeqxbrg'+sci, fl_inter='no',
                  combine='median', sepslits='yes')

os.system('ds9 &')
iraf.sleep(5.0)
for sci in iraf.type(ic.lst_sci, Stdout=1):
    sci = sci.strip()
    iraf.display('txeqxbrg'+sci+'[sci,1]', 1)
    iraf.display('stxeqxbrg'+sci+'[sci,1]', 2)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))