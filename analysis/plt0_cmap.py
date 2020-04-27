#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 9 12:42:05 2020

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy import wcs
import imgscl
# from scipy.misc import bytescale
from PIL import Image
import pandas as pd

# kpc20 = 20.0/4.370/0.03
ifu_h = 5.0 / 0.05
ifu_w = 7.0 / 0.05
# seeing = 0.7 / 0.05


# ----- Directories ----- #
dir_fig = '/data/jlee/DATA/Gemini/Programs/GN-2019A-Q-215/analysis/diagram/'
diT = '/data/jlee/DATA/HLA/McPartland+16/'
diI = '/data/jlee/DATA/HLA/McPartland+16/MACS1752/Images/'
diG = '/data/jlee/DATA/Gemini/Programs/GN-2019A-Q-215/redux4_700/'


# ----- Reading FITS file ----- #
hdr1 = fits.getheader(diG+'cstxeqxbrgN20190611S0257_3D.fits', ext=0)
gra = hdr1['RA'] ; gdec = hdr1['DEC']


# ----- Reading FITS images and creating RGB data ----- #
#img435 = fits.getdata(diI+'alw435NE_drc.fits')
img606 = fits.getdata(diI+'calw606NE_msk.fits')
img814 = fits.getdata(diI+'calw814NE_msk.fits')

cimg = np.zeros((img814.shape[0], img814.shape[1], 3), dtype=float)
cimg[:,:,0] = imgscl.linear(0.5*img814, scale_min=-0.02, scale_max=0.075)   # R
cimg[:,:,1] = imgscl.linear(0.5*(img606+img814), scale_min=-0.02, scale_max=0.15)   # G
cimg[:,:,2] = imgscl.linear(0.5*img606, scale_min=-0.02, scale_max=0.075)   # B


# ----- WCS to XY ----- #
w = wcs.WCS(diI+'calw606NE_msk.fits')
px, py = w.wcs_world2pix(gra, gdec, 1)


# ----- Figure 1 ----- #
fig1 = plt.figure(1, figsize=(6,6))
ax1 = fig1.add_subplot(1,1,1)
ax1.set_position([0.0,0.0,1.0,1.0])
ax1.tick_params(labelleft=False)
ax1.tick_params(labelbottom=False)
plt.tick_params(width=0.0, length=0.0)
# -------------------- #

# Color map
rth = 125.0
plt.imshow(cimg[int(np.round(py)-1-rth):int(np.round(py)-1+rth),
                int(np.round(px)-1-rth):int(np.round(px)-1+rth),:],
                aspect='equal', origin='lower')

# IFU FOV
ang = (85.0 - 54.9894) * (np.pi/180.0) ; x_off = 0.0 ; y_off = 0.0
rot = np.array([[np.cos(ang),-np.sin(ang)],[np.sin(ang),np.cos(ang)]])

spoint = np.array([[-0.5*ifu_w,+0.5*ifu_h],[+0.5*ifu_w,+0.5*ifu_h],
                   [+0.5*ifu_w,-0.5*ifu_h],[-0.5*ifu_w,-0.5*ifu_h]])
epoint = np.array([[+0.5*ifu_w,+0.5*ifu_h],[+0.5*ifu_w,-0.5*ifu_h],
                   [-0.5*ifu_w,-0.5*ifu_h],[-0.5*ifu_w,+0.5*ifu_h]])

for i in np.arange(4):
    srot = np.dot(rot,spoint[i]) ; erot = np.dot(rot,epoint[i])
    ax1.arrow(srot[0]+rth+x_off, srot[1]+rth+y_off, erot[0]-srot[0], erot[1]-srot[1],
              width=2.0, head_width=0., head_length=0., fc='cyan', ec='cyan', alpha=0.7)

# # Seeing circle
# c0 = plt.Circle((225.0,35.0), radius=0.5*seeing, color='w', linewidth=2.5,
#                 linestyle='-', fill=False, alpha=0.9)
# ax1.add_artist(c0)

# ax1.text(215.5, 12.5, '0.7"', ha='left', va='bottom',
#          fontsize=17.5, fontweight='bold', color='w')

# Scale bar
kpc10 = 10.0/4.967/0.05
ax1.arrow(25.0, 35.0, kpc10, 0., width=2.5, head_width=0., head_length=0.,
          fc='yellow', ec='yellow', alpha=0.9)
ax1.text(25.5, 20.0, '10 kpc', fontsize=20.0, fontweight='bold', color='yellow')


# The orientations
x0 = 205.0 ; y0 = 205.0
L = 25.0 ; theta0 = -54.9894*(np.pi/180.0)
ax1.arrow(x0, y0-1.0, -L*np.sin(theta0), L*np.cos(theta0), width=2.5,
          head_width=8.0, head_length=8.0, fc='yellow', ec='yellow', alpha=0.9)
ax1.arrow(x0+1.0, y0, -L*np.cos(theta0), -L*np.sin(theta0), width=2.5,
          head_width=8.0, head_length=8.0, fc='yellow', ec='yellow', alpha=0.9)
ax1.text(175.0, 232.5, 'E', fontsize=20.0, fontweight='bold', color='yellow')
ax1.text(235.0, 222.5, 'N', fontsize=20.0, fontweight='bold', color='yellow')

# # Target name
# ax1.text(15.0, 225.0, 'MACSJ1752-JFG2', fontsize=24.0, ha='left', va='bottom',
#          fontweight='bold', color='w')
# ax1.text(32.5, 205.0, '(z=0.353)', fontsize=23.0, ha='left', va='bottom',
#          fontweight='bold', color='w')

# Saving the figure
plt.savefig(dir_fig+'cmap.pdf')
plt.savefig(dir_fig+'cmap.png', dpi=300)
plt.close()


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))