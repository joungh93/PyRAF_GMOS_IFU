#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 11:18:12 2020

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import glob, os
import init_cfg as ic
from astropy.io import fits
import copy
from astropy.convolution import convolve
from astropy.convolution import Gaussian1DKernel
import tqdm


# ----- Basic parameters ----- #
redshift = 0.353
dist_lum = 1875.5e+6  # pc


# ----- Reading the cube ----- #
fin_cb = 'fcube_3D.fits'

hd0 = fits.getheader(fin_cb, ext=0)
d_sci, h_sci = fits.getdata(fin_cb, ext=1, header=True)
d_var, h_var = fits.getdata(fin_cb, ext=2, header=True)

wav = np.linspace(start=h_sci['CRVAL3'],
                  stop=h_sci['CRVAL3']+(h_sci['NAXIS3']-1)*h_sci['CD3_3'],
                  num=h_sci['NAXIS3'], endpoint=True)

wav_rest = wav/(1.0+redshift)
# d_sci *= (1.0+redshift)
# d_var *= (1.0+redshift)**2.0

# d_snr = d_sci / np.sqrt(d_var)


# ----- Smoothing spectra ----- #

# Wavelength masking range (rest-frame)
wav_msk = np.array([[3720, 3740],  # [OII]3727
#                     [4110, 4130],  # sky
                    [4850, 4875],  # H beta
                    [4950, 4965],  # [OIII]4959
                    [5000, 5015],  # [OIII]5007
                    [6360, 6420],  # sky
                    [6460, 6510],  # sky
                    [6515, 6537],  # sky
                    [6540, 6574],  # [NII] + H alpha
                    [6578, 6592],  # [NII]
                    [6710, 6740]   # [SII]
                   ])

# Copying cube data & creating kernels
d_sci2 = copy.deepcopy(d_sci)
g1 = Gaussian1DKernel(stddev = 40)
g2 = Gaussian1DKernel(stddev = 4)

for x in tqdm.trange(d_sci2.shape[2]):
    for y in np.arange(d_sci2.shape[1]):
        for i in np.arange(wav_msk.shape[0]):
            spx_l = np.abs(wav_rest-wav_msk[i,0]).argmin()
            spx_r = np.abs(wav_rest-wav_msk[i,1]).argmin()
            d_sci2[spx_l:spx_r+1, y, x] = np.nan
        
        d_filt1 = convolve(d_sci2[:, y, x], g1)
        
        for i in np.arange(wav_msk.shape[0]):
            spx_l = np.abs(wav_rest-wav_msk[i,0]).argmin()
            spx_r = np.abs(wav_rest-wav_msk[i,1]).argmin()
            fill_value = d_filt1[spx_l:spx_r+1]
            d_sci2[spx_l:spx_r+1, y, x] = fill_value       
        
        d_filt2 = convolve(d_sci2[:, y, x], g2)
        
        d_sci2[:,y,x] = d_sci[:,y,x] - d_filt2


# ----- Creating the final cube for analysis ----- #
fhd0 = fits.PrimaryHDU()
fhd1 = fits.ImageHDU()
fhd2 = fits.ImageHDU()

fhd0.header = hd0

fhd1.data = d_sci2
fhd1.header = h_sci

fhd2.data = d_var
fhd2.header = h_var

fcb_hdu = fits.HDUList([fhd0, fhd1, fhd2])
fcb_hdu.writeto('sfcube_3D.fits', overwrite=True)


# Printing the running time
print('\n')
print('--- %s seconds ---' %(time.time()-start_time))
