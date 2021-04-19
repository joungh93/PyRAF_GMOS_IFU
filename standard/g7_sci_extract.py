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
from astropy.io import fits


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

iraf.unlearn('gqecorr')
iraf.unlearn('gfextract')
iraf.unlearn('gfdisplay')


# ---------- QE correction and extraction of the science ---------- #
flat = np.loadtxt(ic.lst_flat, dtype=str)
if (flat.size > 1):
    raise ValueError("Please check if there is only one flat image for the standard star.")
flat0 = flat.item(0)

arc = np.loadtxt(ic.lst_arc, dtype=str)
if (arc.size > 1):
    raise ValueError("Please check if there is only one arc image for the standard star.")
arc0 = arc.item(0)

std = np.loadtxt(ic.lst_std, dtype=str)
if (std.size > 1):
    raise ValueError("Please check if there is only one image for the standard star.")
std0 = std.item(0)

response = flat0+'_resp'
ref_flat0 = 'eqbrg'+flat0

iraf.imdelete('qxbrg@'+ic.lst_std)
iraf.imdelete('eqxbrg@'+ic.lst_std)
iraf.gqecorr('xbrg'+std0, refimage='erg'+arc0, fl_correct='yes',
             fl_vardq='yes', verbose='yes')
iraf.gfextract('qxbrg'+std0, response=response, recenter='no',
               trace='no', reference=ref_flat0, weights='none',
               fl_vardq='yes', line=pk_line, exslits=ic.eslit)

if (ic.nslit == 1):
    vkw = '1'
if (ic.nslit == 2):
    vkw = '*'
os.system('ds9 &')
iraf.sleep(5.0)
iraf.gfdisplay('eqxbrg'+std0, 1, version=vkw)


# ----- Displaying the extracted data for checking bad columns ----- #
# fits.open('eqxbrg'+std0+'.fits').info()

dt1, hd1 = fits.getdata('eqxbrg'+std0+'.fits', ext=2, header=True)
z1l, z1u = np.percentile(dt1, [15, 85])
# 2-slit mode
if (ic.nslit == 2):
	dt2, hd2 = fits.getdata('eqxbrg'+std0+'.fits', ext=5, header=True)
	z2l, z2u = np.percentile(dt2, [15, 85])

if (ic.nslit == 1):
	z1, z2 = z1l, z1u
	ds9_frm = "ds9 eqxbrg"+std0+".fits[2] -multiframe"
	ds9_loc = " -scale lock yes -frame lock image"
	ds9_scl = " -scale limits {0:.2f} {1:.2f}".format(z1, z2)
if (ic.nslit == 2):
	z1, z2 = 0.5*(z1l+z2l), 0.5*(z1u+z2u)
	ds9_frm = "ds9 eqxbrg"+std0+".fits[2] eqxbrg"+std0+".fits[5] -multiframe"
	ds9_loc = " -scale lock yes -frame lock image"
	ds9_scl = " -scale limits {0:.2f} {1:.2f}".format(z1, z2)

os.system(ds9_frm + ds9_loc + ds9_scl)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
