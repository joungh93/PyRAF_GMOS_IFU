#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 16:53:20 2019

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import glob, os
from astropy.io import fits
import init_cfg as ic


# ----- Directories ----- #
if (glob.glob(ic.dir_fig+"cubes/") == []):
    os.system("mkdir "+ic.dir_fig+"cubes/")


# ----- Wavelength cutout ranges ----- #
wav_range_obs = (1+ic.redshift) * ic.wav_range_res
f = open("count_max.txt", "w")
for i in np.arange(len(ic.cube_list)):
    hd0 = fits.getheader(ic.cube_list[i], ext=0)
    d_sci, h_sci = fits.getdata(ic.cube_list[i], ext=1, header=True)
    # d_var, h_var = fits.getdata(ic.cube_list[i], ext=2, header=True)

    d_sci[d_sci < 0.0] = 0.0

    wav = np.linspace(start=h_sci['CRVAL3']+(1-h_sci['CRPIX3'])*h_sci['CD3_3'],
                      stop=h_sci['CRVAL3']+(h_sci['NAXIS3']-h_sci['CRPIX3'])*h_sci['CD3_3'],
                      num=h_sci['NAXIS3'], endpoint=True)
    check_Ha = [np.abs(wav-wav_range_obs[0]).argmin(),
                np.abs(wav-wav_range_obs[1]).argmin()+1]

    # ----- Figure 1 : Histogram ----- #
    fig = plt.figure(1, figsize=(9,6))
    ax = plt.subplot(1,1,1)
    ax.set_position([0.12,0.15,0.85,0.82])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xticks([1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0])
    ax.set_xticklabels([r'$10^{-7}$', r'$10^{-6}$', r'$10^{-5}$', r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$'], fontsize=16.0)
    ax.set_yticks([1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0])
    ax.set_yticklabels([r'$10^{0}$', r'$10^{1}$', r'$10^{2}$', r'$10^{3}$', r'$10^{4}$', r'$10^{5}$'], fontsize=16.0)
    ax.set_xlabel('Pixel values', fontsize=16.0)
    ax.set_ylabel(r'$N$', fontsize=16.0)
    ax.set_xlim([1.0e-7,2.0])
    ax.set_ylim([9.0e-1, 5.0e+5])
    ax.tick_params(width=2.0, length=12.0)
    plt.minorticks_on()
    ax.tick_params(width=2.0, length=8.0, which='minor')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2.0)
    # -------------------------------- #

    # count, base = np.histogram(d_sci[:, 2:46, 33:65], bins=np.logspace(-7.0, 0.0, 70))
    ax.hist(d_sci[:, ic.check_y[0]:ic.check_y[1], ic.check_x[0]:ic.check_x[1]].ravel(),
            bins=np.logspace(-7.0, 0.0, 70), color='dodgerblue', alpha=0.8, label="Total wavelength")
    ax.hist(d_sci[check_Ha[0]:check_Ha[1], ic.check_y[0]:ic.check_y[1], ic.check_x[0]:ic.check_x[1]].ravel(),
            bins=np.logspace(-7.0, 0.0, 70), color='magenta', alpha=0.8, label=r"${\rm H\alpha}$ region")
    hist, bin_edges = np.histogram(d_sci[check_Ha[0]:check_Ha[1], ic.check_y[0]:ic.check_y[1], ic.check_x[0]:ic.check_x[1]].ravel(),
                                   bins=np.logspace(-7.0, 0.0, 70))
    hist = hist.astype('float')
    hist[hist == 0] = 0.9
    log_hist = np.log10(hist)
    dlog_hist = np.diff(log_hist)
    idx_max = np.abs(dlog_hist).argmax()
    plt_max = 10**(0.5*(np.log10(bin_edges[idx_max]) + np.log10(bin_edges[idx_max+1])))
    count_max = bin_edges[idx_max+1]
    # count_max = np.max(d_sci[check_Ha[0]:check_Ha[1], ic.check_y[0]:ic.check_y[1], ic.check_x[0]:ic.check_x[1]])
    ax.axvline(plt_max, 0, 1, color='k', linestyle='--', linewidth=2.0, alpha=0.8)
    
    hds, lbs = ax.get_legend_handles_labels()
    ax.legend(hds, lbs, fontsize=12.0, loc=(0.68, 0.82),
              handlelength=2.5, frameon=True, borderpad=0.8,
              framealpha=0.8, edgecolor='gray')

    f.write(ic.cube_name[i]+f"\t{count_max:.4e}\n")
    d_sci[d_sci > count_max] = 0.0

    plt.savefig(ic.dir_fig+'cubes/'+'hist-'+ic.cube_name[i]+'.png', dpi=300)
    plt.close()


    # ----- Figure 2 : wavelength cutout images ----- #
    fig = plt.figure(2, figsize=(10.2, 7.35))
    ax = plt.subplot(1,1,1)
    ax.set_position([0.,0.,1.,1.])
    ax.tick_params(labelleft=False)
    ax.tick_params(labelbottom=False)
    ax.tick_params(width=0.0, length=0.0)
    # ----------------------------------------------- #
    d_white = d_sci[check_Ha[0]:check_Ha[1],:,:].sum(axis=0)
    v_low, v_high = np.percentile(d_white, [5.0, 95.0])
    ax.imshow(d_white, cmap='gray', vmin=v_low, vmax=v_high, origin='lower')
    p0 = Rectangle((ic.check_x[0], ic.check_y[0]), width = ic.check_x[1]-ic.check_x[0],
                   height = ic.check_y[1]-ic.check_y[0], angle = 0.0,
                   ls = "-", lw = 4, ec = "orange", fill = False, alpha = 0.7)
    ax.add_patch(p0)

    plt.savefig(ic.dir_fig+'cubes/'+'sum-'+ic.cube_name[i]+'_'+ \
                str(round(wav_range_obs[0]))+'-'+str(round(wav_range_obs[1]))+'A.png', dpi=300)
    plt.close()

    # Saving H alpha FITS images
    fits.writeto(ic.dir_cmb+'Ha_sum-'+ic.cube_name[i]+'.fits', d_white, h_sci, overwrite=True)
f.close()


# Printing the running time
print('--- %.4f seconds ---' %(time.time()-start_time))
