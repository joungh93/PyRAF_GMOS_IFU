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


# ----- Importing IRAF from the root directory ----- #
current_dir = os.getcwd()
os.chdir(ic.dir_iraf)

from pyraf import iraf
from pyraf.iraf import gemini, gmos

os.chdir(current_dir)
iraf.chdir(current_dir)

iraf.unlearn('gfapsum')
iraf.unlearn('gsstandard')


# ---------- Sensitivity function (after reducing standard star) ---------- #
std = np.loadtxt(ic.lst_std, dtype=str)
if (std.size > 1):
    raise ValueError("Please check if there is only one image for the standard star.")
std0 = std.item(0)

# Sum the fibers
iraf.imdelete('astxeqxbrg@'+ic.lst_std)

iraf.gfapsum('stxeqxbrg'+std0, combine='sum', fl_inter='no')

for std in iraf.type(ic.lst_std, Stdout=1):
    std = std.strip()
    fits.open('astxeqxbrg'+std+'.fits').info()
    iraf.splot('astxeqxbrg'+std+'[sci,1]')


# Call gsstandard
outflux = ic.root_name+'std'
sensfunc = ic.root_name+'sens'

iraf.delete(outflux, verify='no')    # Not .fits file
iraf.imdelete(sensfunc, verify='no')    # .fits file

for std in iraf.type(ic.lst_std, Stdout=1):
    std = std.strip()
    iraf.gsstandard('astxeqxbrg'+std, outflux, sensfunc,
    	            starname = ic.starname, observatory = ic.obs_site,
    	            caldir = ic.stardir, extinction = ic.extinction,
    	            fl_inter = 'yes', function = 'chebyshev', order=11)

# Store the solution
iraf.copy(sensfunc+'.fits', ic.caldir)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
