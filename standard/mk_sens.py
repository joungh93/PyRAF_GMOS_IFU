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
from matplotlib import pyplot as plt


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
iraf.splot('astxeqxbrg'+std0+'[sci,1]')


# Plot the sum spectra
spec_sum, hd = fits.getdata('astxeqxbrg'+std0+'.fits', ext=1, header=True)
w = hd['CRVAL1'] + np.arange(len(spec_sum))*hd['CD1_1']

fig, ax = plt.subplots(1, 1, figsize=(7,5))
plt.subplots_adjust(left=0.16, right=0.96, bottom=0.13, top=0.91)
plt.suptitle("Sum of spectra", x=0.55, ha='center', y=0.97, va='top',
             fontsize=16.0, fontweight='bold')
ax.set_xlabel(r"Wavelength [${\rm \AA}$]", fontsize=12.0)
ax.set_ylabel("Counts", fontsize=12.0)
# ax.set_yscale('log')
ax.tick_params(axis='both', labelsize=12.0)
ax.plot(w, spec_sum, color='C0', linewidth=1.5, alpha=0.8)
ax.text(0.96, 0.95, ic.starname, fontsize=14.0, fontweight='bold',
         color='k', ha='right', va='top', transform=ax.transAxes)
plt.savefig(ic.caldir+ic.starname+"_specsum.png", dpi=300)
plt.close()


# Call gsstandard
outflux = ic.root_name+'std'
sensfunc = ic.root_name+'sens'

iraf.delete(outflux, verify='no')    # Not .fits file
iraf.imdelete(sensfunc, verify='no')    # .fits file

iraf.gsstandard('astxeqxbrg'+std0, outflux, sensfunc,
                starname = ic.starname, observatory = ic.obs_site,
                caldir = ic.stardir, extinction = ic.extinction,
                fl_inter = 'yes', function = 'chebyshev', order=11)

# Store the solution
os.system("cp -rpv "+sensfunc+".fits "+ic.caldir)


# Printing the running time
print('--- %.3f seconds ---' %(time.time()-start_time))
