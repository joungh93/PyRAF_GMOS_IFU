#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 9 14:45:23 2020

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import glob, os, copy
import init_cfg as ic
from astropy.io import fits
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import tqdm
from scipy.special import erf
import pandas as pd
from scipy import ndimage
from astropy.cosmology import FlatLambdaCDM
import warnings
warnings.filterwarnings("ignore")


# ----- Directory ----- #
dir_fig = '/data/jlee/DATA/Gemini/Programs/GN-2019A-Q-215/analysis/diagram/linefits/'
diG = '/data/jlee/DATA/Gemini/Programs/GN-2019A-Q-215/redux4_700/'
dir_lines = '/data/jlee/DATA/Gemini/Programs/GN-2019A-Q-215/analysis/lines2/'
glob_lines = glob.glob(dir_lines+'*')

emi_lines = []
for i in np.arange(len(glob_lines)):
    emi_lines.append(glob_lines[i].split(dir_lines)[1])
emi_lines = sorted(emi_lines)
emi_lines.remove('check')


for i in np.arange(len(emi_lines)):
    exec('dir_'+emi_lines[i]+' = "'+dir_lines+emi_lines[i]+'/"')

# emi_lines
# ['Halpha',
#  'Hbeta',
#  'NII6548',
#  'NII6584',
#  'OIII4959',
#  'OIII5007',
#  'SII6717',
#  'SII6731']

name_elines = [r"${\rm H\alpha}$",
               r"${\rm H\beta}$",
               r"${\rm [NII]\lambda6548}$",
               r"${\rm [NII]\lambda6584}$",
               r"${\rm [OIII]\lambda4959}$",
               r"${\rm [OIII]\lambda5007}$",
               r"${\rm [SII]\lambda6717}$",
               r"${\rm [SII]\lambda6731}$"]


# ----- Basic parameters ----- #
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
redshift = 0.3527
dist_lum = cosmo.luminosity_distance(redshift).value*1.0e+6    # pc
c = 2.99792e+5    # km/s
ang_scale = cosmo.kpc_proper_per_arcmin(redshift).value / 60.    # kpc/arcsec
pixel_scale = 0.1    # arcsec/pix

# Angstrom (SDSS)
wav_lines = [6564.61, 4862.68, 6549.86, 6585.27,
             4960.295, 5008.240, 6718.29, 6732.67]
wav_lines = np.array(wav_lines)

# Spectral resolution & sigma0
par, e_par = np.loadtxt('relation_wav_R.txt').T
Rmax = par[0] + par[1]*wav_lines[0]*(1+redshift)
e_Rmax = np.sqrt(e_par[0]**2.0 + (e_par[1]*wav_lines[0]*(1+redshift))**2.0)
lsig0 = wav_lines[0] / (2.0*np.sqrt(2.0*np.log(2.0))*Rmax)
e_lsig0 = wav_lines[0]*e_Rmax / (2.0*np.sqrt(2.0*np.log(2.0))*Rmax*Rmax)
lsig_llim = lsig0# - 1.0*e_lsig0
vsig0 = c / (2.0*np.sqrt(2.0*np.log(2.0))*Rmax)

for i in np.arange(len(emi_lines)):
    exec('wav_'+emi_lines[i]+' = '+str(wav_lines[i]))


# ----- Reading 2D data ----- #
for l in np.arange(len(emi_lines)):
    exec("dat_dir = dir_"+emi_lines[l])
    dat2D_list = sorted(glob.glob(dat_dir+'*_2D.fits'))
    for i in np.arange(len(dat2D_list)):
        dat_2D = emi_lines[l]+'_'+dat2D_list[i].split('/')[-1].split('.fits')[0]
        exec(dat_2D+" = fits.getdata('"+dat2D_list[i]+"', ext=0)")

for l in np.arange(len(emi_lines)):
    exec("flux = "+emi_lines[l]+"_flux_2D")
    exec("snrpix = "+emi_lines[l]+"_snrpix_2D")
    fluxerr = flux / snrpix
    exec("e_"+emi_lines[l]+"_flux_2D = fluxerr")

# Reading central RA & Dec
hdr1 = fits.getheader(diG+'cstxeqxbrgN20190611S0257_3D.fits', ext=0)
gra, gdec, gpa = hdr1['RA'], hdr1['DEC'], hdr1['PA']


# ----- Loading Voronoi binned data ----- #
dir_vbin = 'vorbin/'
vb = np.load(dir_vbin+'vorbin_array.npz')
# wav, sci, var
data_vbin = fits.getdata(dir_vbin+'vbin.fits').astype('int')
nvbin = np.unique(data_vbin).size-1

N_area = np.array([], dtype='int')
for ibin in np.arange(nvbin):
    nvbin_region = (data_vbin == ibin)
    N_area = np.append(N_area, np.sum(nvbin_region))


# ----- Creating H alpha contour data ----- #
X_coord = 0.1*(np.arange(Halpha_flux_2D.shape[1], step=1)-Halpha_flux_2D.shape[1]/2)
Y_coord = 0.1*(np.arange(Halpha_flux_2D.shape[0], step=1)-Halpha_flux_2D.shape[0]/2)

flx_Data = copy.deepcopy(Halpha_flux_2D)
snr_cnd = ((Halpha_snr_2D < 3.0) | (Halpha_snrpix_2D < 5.0))
sig_cnd = (Halpha_sigma_2D < lsig_llim)
rchi25, rchi50, rchi75 = np.percentile(Halpha_rchisq_2D[Halpha_rchisq_2D > 0.],
                                       [25.0, 50.0, 75.0])
rchisq_cnd = (Halpha_rchisq_2D > 50.)
flx25, flx50, flx75 = np.percentile(Halpha_flux_2D[Halpha_flux_2D > 0.],
                                    [25.0, 50.0, 75.0])
flx_cnd = (Halpha_flux_2D < flx25)
# zero_cnd = (snr_cnd | rchisq_cnd)
zero_cnd = (snr_cnd | rchisq_cnd)# | flx_cnd)

flx_Data[zero_cnd] = 0.
sflx = ndimage.gaussian_filter(flx_Data, sigma=(1.2,1.2), order=0)

sig = np.std(sflx)
lvs = [0.5*sig, 1.*sig, 2.*sig, 3.*sig, 5.*sig, 7.*sig]
lws = tuple(np.repeat(2.5, len(lvs)-1))
cs = tuple(['gray']*(len(lvs)-1))


# ----- START: Flux maps ----- #
for l in np.arange(len(emi_lines)):
    exec("snr_cnd = ((Halpha_snr_2D < 3.0) | (Halpha_snrpix_2D < 5.0) | ("+emi_lines[l]+"_snr_2D < 3.0))")
    sig_cnd = (Halpha_sigma_2D < lsig_llim)
    rchisq_cnd = (Halpha_rchisq_2D > 50.)
    flx_cnd = (Halpha_flux_2D < flx25)
    zero_cnd = (snr_cnd | rchisq_cnd)

    exec(emi_lines[l]+"_flux_2D[zero_cnd] = 0.")
    exec("plt_Data = "+emi_lines[l]+"_flux_2D")

    fig, ax = plt.subplots(1, 1, figsize=(8,5))
    plt.suptitle(name_elines[l]+" flux map",
                 x=0.5, ha='center', y=0.96, va='top',
                 fontsize=20.0)
    ax.set_xlim([-3.4, 3.4])
    ax.set_ylim([-2.45, 2.45])
    ax.set_xticks([-3,-2,-1,0,1,2,3])
    ax.set_yticks([-2,-1,0,1,2])
    ax.set_xticklabels([r'$-3$',r'$-2$',r'$-1$',0,1,2,3], fontsize=15.0)
    ax.set_yticklabels([r'$-2$',r'$-1$',0,1,2], fontsize=15.0)
    ax.set_xlabel('arcsec', fontsize=15.0) 
    ax.set_ylabel('arcsec', fontsize=15.0)
    ax.tick_params(width=1.0, length=5.0)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.0)

    v_low, v_high = np.percentile(plt_Data[plt_Data > 0.], [1.0, 99.0])
    im = ax.imshow(plt_Data, cmap='gray_r', vmin=0., vmax=v_high,
                   aspect='equal', extent=[-3.4,3.4,-2.45,2.45])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(im, cax=cax)
    cb.set_label(r'Flux [${\rm 10^{-15}~erg~s^{-1}~cm^{-2}~\AA^{-1}}$]',
                 size=12.0, labelpad=15.0)
    cb.ax.tick_params(labelsize=12.0)

    ax.contour(X_coord, Y_coord[::-1], sflx, levels=lvs, linewidths=lws, colors=cs, alpha=0.6)
    p0, = ax.plot(-100.0, -100.0, '-', linewidth=2.5, color='gray', alpha=0.6,
                  label=r"H${\rm \alpha}$ flux contour")
    ax.legend(handles=[p0], fontsize=13.0, loc='lower left',
              handlelength=2.5, frameon=True, borderpad=0.8,
              framealpha=0.8, edgecolor='gray')

    # The orientations
    x0 = -2.75 ; y0 = 1.25
    L = 0.6 ; theta0 = gpa*(np.pi/180.0)
    ax.arrow(x0-0.025, y0, L*np.sin(theta0), L*np.cos(theta0), width=0.06,
             head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
    ax.arrow(x0, y0-0.025, -L*np.cos(theta0), L*np.sin(theta0), width=0.06,
             head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
    ax.text(-2.95, 2.10, 'E', fontsize=15.0, fontweight='bold', color='blueviolet')
    ax.text(-1.90, 1.25, 'N', fontsize=15.0, fontweight='bold', color='blueviolet')

    # Scale bar
    kpc5 = 5.0 / ang_scale
    ax.arrow(2.0, -1.85, kpc5, 0., width=0.07, head_width=0., head_length=0.,
              fc='blueviolet', ec='blueviolet', alpha=0.9)
    ax.text(2.1, -2.2, '5 kpc', fontsize=15.0, fontweight='bold', color='blueviolet')

    plt.savefig(dir_fig+'Map_flux_'+emi_lines[l]+'.pdf')
    plt.savefig(dir_fig+'Map_flux_'+emi_lines[l]+'.png', dpi=300)
    plt.close()
# ----- END: Flux maps ----- #


# ----- START: S/N maps ----- #
for l in np.arange(len(emi_lines)):
    # exec("snr_pix = copy.deepcopy("+emi_lines[l]+"_snr_2D)")
    # for ibin in np.arange(nvbin):
    #     nvbin_region = (data_vbin == ibin)
    #     exec("snr_pix[nvbin_region] = "+emi_lines[l]+"_snr_2D[nvbin_region]/np.sqrt(N_area[ibin])")
    
    # exec("snr_cnd = ((Halpha_snr_2D < 3.0) | ("+emi_lines[l]+"_snr_2D < 3.0))")
    # sig_cnd = (Halpha_sigma_2D < lsig_llim)
    # rchisq_cnd = (Halpha_rchisq_2D > 50.)
    # flx_cnd = (Halpha_flux_2D < flx25)
    # zero_cnd = (snr_cnd | rchisq_cnd)

    # snr_pix[zero_cnd] = 0.
    # plt_Data = snr_pix
    exec("plt_Data = "+emi_lines[l]+"_snrpix_2D")

    fig, ax = plt.subplots(1, 1, figsize=(8,5))
    plt.suptitle(name_elines[l]+" S/N map (per pixel)",
                 x=0.5, ha='center', y=0.96, va='top',
                 fontsize=20.0)
    ax.set_xlim([-3.4, 3.4])
    ax.set_ylim([-2.45, 2.45])
    ax.set_xticks([-3,-2,-1,0,1,2,3])
    ax.set_yticks([-2,-1,0,1,2])
    ax.set_xticklabels([r'$-3$',r'$-2$',r'$-1$',0,1,2,3], fontsize=15.0)
    ax.set_yticklabels([r'$-2$',r'$-1$',0,1,2], fontsize=15.0)
    ax.set_xlabel('arcsec', fontsize=15.0) 
    ax.set_ylabel('arcsec', fontsize=15.0)
    ax.tick_params(width=1.0, length=5.0)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.0)

    v_low, v_high = np.percentile(plt_Data[plt_Data > 0.], [1.0, 99.0])
    im = ax.imshow(plt_Data, cmap='gray_r', vmin=0., vmax=v_high,
                   aspect='equal', extent=[-3.4,3.4,-2.45,2.45])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(im, cax=cax)
    cb.set_label('Signal-to-noise ratio per pixel',
                 size=12.0, labelpad=15.0)
    cb.ax.tick_params(labelsize=12.0)\

    ax.contour(X_coord, Y_coord[::-1], sflx, levels=lvs, linewidths=lws, colors=cs, alpha=0.6)
    p0, = ax.plot(-100.0, -100.0, '-', linewidth=2.5, color='gray', alpha=0.6,
                  label=r"H${\rm \alpha}$ flux contour")
    ax.legend(handles=[p0], fontsize=13.0, loc='lower left',
              handlelength=2.5, frameon=True, borderpad=0.8,
              framealpha=0.8, edgecolor='gray')

    # The orientations
    x0 = -2.75 ; y0 = 1.25
    L = 0.6 ; theta0 = gpa*(np.pi/180.0)
    ax.arrow(x0-0.025, y0, L*np.sin(theta0), L*np.cos(theta0), width=0.06,
             head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
    ax.arrow(x0, y0-0.025, -L*np.cos(theta0), L*np.sin(theta0), width=0.06,
             head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
    ax.text(-2.95, 2.10, 'E', fontsize=15.0, fontweight='bold', color='blueviolet')
    ax.text(-1.90, 1.25, 'N', fontsize=15.0, fontweight='bold', color='blueviolet')

    # Scale bar
    kpc5 = 5.0 / ang_scale
    ax.arrow(2.0, -1.85, kpc5, 0., width=0.07, head_width=0., head_length=0.,
              fc='blueviolet', ec='blueviolet', alpha=0.9)
    ax.text(2.1, -2.2, '5 kpc', fontsize=15.0, fontweight='bold', color='blueviolet')

    plt.savefig(dir_fig+'Map_snr_'+emi_lines[l]+'.pdf')
    plt.savefig(dir_fig+'Map_snr_'+emi_lines[l]+'.png', dpi=300)
    plt.close()
# ----- END: S/N maps ----- #


# ----- START: Radial velocity distribution (H alpha) map ----- #
ymax_idx, xmax_idx = np.unravel_index(Halpha_flux_2D.argmax(), Halpha_flux_2D.shape)
Halpha_mu0 = Halpha_mu_2D[ymax_idx, xmax_idx]

snr_cnd = ((Halpha_snr_2D < 3.0) | (Halpha_snrpix_2D < 5.0))
sig_cnd = (Halpha_sigma_2D < lsig_llim)
rchisq_cnd = (Halpha_rchisq_2D > 50.)
flx_cnd = (Halpha_flux_2D < flx25)
zero_cnd = (snr_cnd | rchisq_cnd)

plt_Data = (c*(Halpha_mu_2D-Halpha_mu0)/Halpha_mu0)
plt_Data[Halpha_mu_2D == 0] = np.nan
plt_Data[zero_cnd] = np.nan
Rvd = plt_Data

fig, ax = plt.subplots(1, 1, figsize=(8,5))
plt.suptitle("Radial velocity map",
             x=0.5, ha='center', y=0.96, va='top',
             fontsize=20.0)
ax.set_xlim([-3.4, 3.4])
ax.set_ylim([-2.45, 2.45])
ax.set_xticks([-3,-2,-1,0,1,2,3])
ax.set_yticks([-2,-1,0,1,2])
ax.set_xticklabels([r'$-3$',r'$-2$',r'$-1$',0,1,2,3], fontsize=15.0)
ax.set_yticklabels([r'$-2$',r'$-1$',0,1,2], fontsize=15.0)
ax.set_xlabel('arcsec', fontsize=15.0) 
ax.set_ylabel('arcsec', fontsize=15.0)
ax.tick_params(width=1.0, length=5.0)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
im = ax.imshow(plt_Data, cmap='rainbow',
               vmin=v_low-25.0, vmax=v_high+25.0, 
               aspect='equal', extent=[-3.4,3.4,-2.45,2.45])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(im, cax=cax)
cb.set_label(r'Relative velocity [${\rm km~s^{-1}}$]',
             size=12.0, labelpad=15.0)
cb.ax.tick_params(labelsize=12.0)

ax.contour(X_coord, Y_coord[::-1], sflx, levels=lvs, linewidths=lws, colors=cs, alpha=0.6)
p0, = ax.plot(-100.0, -100.0, '-', linewidth=2.5, color='gray', alpha=0.6,
              label=r"H${\rm \alpha}$ flux contour")
ax.legend(handles=[p0], fontsize=13.0, loc='lower left',
          handlelength=2.5, frameon=True, borderpad=0.8,
          framealpha=0.8, edgecolor='gray')

# The orientations
x0 = -2.75 ; y0 = 1.25
L = 0.6 ; theta0 = gpa*(np.pi/180.0)
ax.arrow(x0-0.025, y0, L*np.sin(theta0), L*np.cos(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.arrow(x0, y0-0.025, -L*np.cos(theta0), L*np.sin(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(-2.95, 2.10, 'E', fontsize=15.0, fontweight='bold', color='blueviolet')
ax.text(-1.90, 1.25, 'N', fontsize=15.0, fontweight='bold', color='blueviolet')

# Scale bar
kpc5 = 5.0 / ang_scale
ax.arrow(2.0, -1.85, kpc5, 0., width=0.07, head_width=0., head_length=0.,
          fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(2.1, -2.2, '5 kpc', fontsize=15.0, fontweight='bold', color='blueviolet')

plt.savefig(dir_fig+'Map_rv_Halpha.pdf')
plt.savefig(dir_fig+'Map_rv_Halpha.png', dpi=300)
plt.close()
# ----- END: Radial velocity distribution (H alpha) map ----- #


# ----- START: Velocity dispersion (H alpha) map ----- #
plt_Data = np.sqrt(Halpha_vsig_2D**2.0 - vsig0**2.0)

snr_cnd = ((Halpha_snr_2D < 3.0) | (Halpha_snrpix_2D < 5.0))
sig_cnd = (Halpha_sigma_2D < lsig_llim)
rchisq_cnd = (Halpha_rchisq_2D > 50.)
flx_cnd = (Halpha_flux_2D < flx25)
zero_cnd = (snr_cnd | rchisq_cnd)

plt_Data[zero_cnd] = np.nan
Vdd = plt_Data

fig, ax = plt.subplots(1, 1, figsize=(8,5))
plt.suptitle("Velocity dispersion map",
             x=0.5, ha='center', y=0.96, va='top',
             fontsize=20.0)
ax.set_xlim([-3.4, 3.4])
ax.set_ylim([-2.45, 2.45])
ax.set_xticks([-3,-2,-1,0,1,2,3])
ax.set_yticks([-2,-1,0,1,2])
ax.set_xticklabels([r'$-3$',r'$-2$',r'$-1$',0,1,2,3], fontsize=15.0)
ax.set_yticklabels([r'$-2$',r'$-1$',0,1,2], fontsize=15.0)
ax.set_xlabel('arcsec', fontsize=15.0) 
ax.set_ylabel('arcsec', fontsize=15.0)
ax.tick_params(width=1.0, length=5.0)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
im = ax.imshow(plt_Data, cmap='rainbow',
               vmin=0.0, vmax=v_high+25.0, 
               aspect='equal', extent=[-3.4,3.4,-2.45,2.45])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(im, cax=cax)
cb.set_label(r'Velocity dispersion [${\rm km~s^{-1}}$]',
             size=12.0, labelpad=15.0)
cb.ax.tick_params(labelsize=12.0)

ax.contour(X_coord, Y_coord[::-1], sflx, levels=lvs, linewidths=lws, colors=cs, alpha=0.6)
p0, = ax.plot(-100.0, -100.0, '-', linewidth=2.5, color='gray', alpha=0.6,
              label=r"H${\rm \alpha}$ flux contour")
ax.legend(handles=[p0], fontsize=13.0, loc='lower left',
          handlelength=2.5, frameon=True, borderpad=0.8,
          framealpha=0.8, edgecolor='gray')

# The orientations
x0 = -2.75 ; y0 = 1.25
L = 0.6 ; theta0 = gpa*(np.pi/180.0)
ax.arrow(x0-0.025, y0, L*np.sin(theta0), L*np.cos(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.arrow(x0, y0-0.025, -L*np.cos(theta0), L*np.sin(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(-2.95, 2.10, 'E', fontsize=15.0, fontweight='bold', color='blueviolet')
ax.text(-1.90, 1.25, 'N', fontsize=15.0, fontweight='bold', color='blueviolet')

# Scale bar
kpc5 = 5.0 / ang_scale
ax.arrow(2.0, -1.85, kpc5, 0., width=0.07, head_width=0., head_length=0.,
          fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(2.1, -2.2, '5 kpc', fontsize=15.0, fontweight='bold', color='blueviolet')

plt.savefig(dir_fig+'Map_vd_Halpha.pdf')
plt.savefig(dir_fig+'Map_vd_Halpha.png', dpi=300)
plt.close()

fwm_vdisp = np.average(plt_Data[np.isnan(plt_Data) == False],
                       weights=flx_Data[np.isnan(plt_Data) == False])

print(f"Flux-weighted mean of velocity dispersion : {fwm_vdisp:.2f} km/s")
# ----- END: Velocity dispersion (H alpha) map ----- #


# ----- START: H alpha / H beta flux ratio map ----- #
plt_Data = Halpha_flux_2D / Hbeta_flux_2D
plt_Data[plt_Data == 0.] = np.nan
plt_Data[np.isinf(plt_Data) == True] = np.nan

snr_cnd = ((Halpha_snr_2D < 3.0) | (Halpha_snrpix_2D < 5.0) | (Hbeta_snr_2D < 3.0))
sig_cnd = (Halpha_sigma_2D < lsig_llim)
rchisq_cnd = (Halpha_rchisq_2D > 50.)
flx_cnd = (Halpha_flux_2D < flx25)
zero_cnd = (snr_cnd | rchisq_cnd)

plt_Data[zero_cnd] = np.nan

fig, ax = plt.subplots(1, 1, figsize=(8,5))
plt.suptitle(r"${\rm H\alpha/H\beta}$ flux ratio map",
             x=0.5, ha='center', y=0.96, va='top',
             fontsize=20.0)
ax.set_xlim([-3.4, 3.4])
ax.set_ylim([-2.45, 2.45])
ax.set_xticks([-3,-2,-1,0,1,2,3])
ax.set_yticks([-2,-1,0,1,2])
ax.set_xticklabels([r'$-3$',r'$-2$',r'$-1$',0,1,2,3], fontsize=15.0)
ax.set_yticklabels([r'$-2$',r'$-1$',0,1,2], fontsize=15.0)
ax.set_xlabel('arcsec', fontsize=15.0) 
ax.set_ylabel('arcsec', fontsize=15.0)
ax.tick_params(width=1.0, length=5.0)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
im = ax.imshow(plt_Data, cmap='rainbow',
               vmin=2.86, vmax=1.1*v_high, 
               aspect='equal', extent=[-3.4,3.4,-2.45,2.45])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(im, cax=cax)

ax.contour(X_coord, Y_coord[::-1], sflx, levels=lvs, linewidths=lws, colors=cs, alpha=0.6)
p0, = ax.plot(-100.0, -100.0, '-', linewidth=2.5, color='gray', alpha=0.6,
              label=r"H${\rm \alpha}$ flux contour")
ax.legend(handles=[p0], fontsize=13.0, loc='lower left',
          handlelength=2.5, frameon=True, borderpad=0.8,
          framealpha=0.8, edgecolor='gray')

# The orientations
x0 = -2.75 ; y0 = 1.25
L = 0.6 ; theta0 = gpa*(np.pi/180.0)
ax.arrow(x0-0.025, y0, L*np.sin(theta0), L*np.cos(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.arrow(x0, y0-0.025, -L*np.cos(theta0), L*np.sin(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(-2.95, 2.10, 'E', fontsize=15.0, fontweight='bold', color='blueviolet')
ax.text(-1.90, 1.25, 'N', fontsize=15.0, fontweight='bold', color='blueviolet')

# Scale bar
kpc5 = 5.0 / ang_scale
ax.arrow(2.0, -1.85, kpc5, 0., width=0.07, head_width=0., head_length=0.,
          fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(2.1, -2.2, '5 kpc', fontsize=15.0, fontweight='bold', color='blueviolet')

plt.savefig(dir_fig+'Line_ratio_Hab.pdf')
plt.savefig(dir_fig+'Line_ratio_Hab.png', dpi=300)
plt.close()

val = (np.isnan(plt_Data) == False)
fwm_Hab = np.average(plt_Data[val], weights=flx_Data[val])

# Error propagation
e_Hab = plt_Data * np.sqrt((e_Halpha_flux_2D/Halpha_flux_2D)**2 + (e_Hbeta_flux_2D/Hbeta_flux_2D)**2)

flx_sum = np.sum(flx_Data[val])
e_flx_sum = np.sqrt(np.sum(e_Halpha_flux_2D[val]**2.))

Aj = (flx_Data[val]*plt_Data[val])
Cj = Aj/flx_sum

e_Aj = Aj * np.sqrt((e_Halpha_flux_2D[val]/flx_Data[val])**2 + (e_Hab[val]/plt_Data[val])**2)
e_Cj = Cj * np.sqrt((e_Aj/Aj)**2. + (e_flx_sum/flx_sum)**2.)

e_fwm_Hab = np.sum(e_Cj)

print(f"Flux-weighted mean of Ha/Hb flux ratio : {fwm_Hab:.3f} +/- {e_fwm_Hab:.3f}")
# ----- END: H alpha / H beta flux ratio map ----- #


# ----- START: SFR (H alpha) map ----- #
# Hab = copy.deepcopy(plt_Data)
Hab = np.ones((plt_Data.shape[0], plt_Data.shape[1]))*fwm_Hab
EBV_gal = 0.026
EBV_int_2D = 1.97*np.log10(Hab/2.86)
k_Halpha = 3.32
k_V = 4.05
A_Halpha = k_Halpha * (EBV_int_2D + EBV_gal)
e_A_Halpha = k_Halpha * 1.97*e_fwm_Hab/(Hab*np.log(10.))
A_V = (k_V/k_Halpha) * A_Halpha

L_Halpha = 1.0e-15 * flx_Data * 10.0**(0.4*A_Halpha) * (4.0*np.pi*(dist_lum*3.086e+18)**2.0)
e_L_Halpha = L_Halpha * np.sqrt((e_Halpha_flux_2D/flx_Data)**2. + (e_A_Halpha/A_Halpha)**2.)
SFR = L_Halpha * 4.6e-42
e_SFR = e_L_Halpha * 4.6e-42
SFR[np.isnan(SFR) == True] = 0.
e_SFR[SFR == 0.] = 0.

plt_Data = SFR / (pixel_scale*ang_scale)**2.
SFRD = plt_Data

fig, ax = plt.subplots(1, 1, figsize=(8,5))
plt.suptitle("SFR map",
             x=0.5, ha='center', y=0.96, va='top',
             fontsize=20.0)
ax.set_xlim([-3.4, 3.4])
ax.set_ylim([-2.45, 2.45])
ax.set_xticks([-3,-2,-1,0,1,2,3])
ax.set_yticks([-2,-1,0,1,2])
ax.set_xticklabels([r'$-3$',r'$-2$',r'$-1$',0,1,2,3], fontsize=15.0)
ax.set_yticklabels([r'$-2$',r'$-1$',0,1,2], fontsize=15.0)
ax.set_xlabel('arcsec', fontsize=15.0) 
ax.set_ylabel('arcsec', fontsize=15.0)
ax.tick_params(width=1.0, length=5.0)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

v_low, v_high = np.percentile(plt_Data, [1.0, 99.0])
im = ax.imshow(plt_Data, cmap='gray_r', vmin=v_low, vmax=v_high,
               aspect='equal', extent=[-3.4,3.4,-2.45,2.45])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(im, cax=cax)
cb.set_label(r'Star formation rate density [$M_{\odot}~{\rm yr^{-1}~kpc^{-2}}$]',
             size=12.0, labelpad=15.0)
cb.ax.tick_params(labelsize=12.0)

ax.contour(X_coord, Y_coord[::-1], sflx, levels=lvs, linewidths=lws, colors=cs, alpha=0.6)
p0, = ax.plot(-100.0, -100.0, '-', linewidth=2.5, color='gray', alpha=0.6,
              label=r"H${\rm \alpha}$ flux contour")
ax.legend(handles=[p0], fontsize=13.0, loc='lower left',
          handlelength=2.5, frameon=True, borderpad=0.8,
          framealpha=0.8, edgecolor='gray')

# The orientations
x0 = -2.75 ; y0 = 1.25
L = 0.6 ; theta0 = gpa*(np.pi/180.0)
ax.arrow(x0-0.025, y0, L*np.sin(theta0), L*np.cos(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.arrow(x0, y0-0.025, -L*np.cos(theta0), L*np.sin(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(-2.95, 2.10, 'E', fontsize=15.0, fontweight='bold', color='blueviolet')
ax.text(-1.90, 1.25, 'N', fontsize=15.0, fontweight='bold', color='blueviolet')

# Scale bar
kpc5 = 5.0 / ang_scale
ax.arrow(2.0, -1.85, kpc5, 0., width=0.07, head_width=0., head_length=0.,
          fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(2.1, -2.2, '5 kpc', fontsize=15.0, fontweight='bold', color='blueviolet')

plt.savefig(dir_fig+'SFR_Halpha.pdf')
plt.savefig(dir_fig+'SFR_Halpha.png', dpi=300)
plt.close()

print(f"SFR sum : {np.sum(SFR):.2f} +/- {np.sqrt(np.sum(e_SFR)):.2f} Mo/yr")
print(f"V-magnitude extinction : {np.mean(A_V):.3f} mag")
# ----- END: SFR (H alpha) map ----- #r


# ----- START: calculating tail SFR ----- #
from photutils.aperture import EllipticalAperture as EAp
from photutils import aperture_photometry as apphot
from reg_saoimage import read_region

def cal_sfr_disk(regfile):
    x0, y0, a, b, pa = read_region(regfile, regtype='ellipse')
    x0, y0 = x0[0]-1, y0[0]-1
    a, b = a[0], b[0]
    pa = pa[0]*np.pi/180.

    ap = EAp(positions=[(x0, y0)], a=a, b=b, theta=pa)
    phot_table = apphot(data=SFR, apertures=ap, error=e_SFR, mask=None)

    return [phot_table['aperture_sum'].data[0], phot_table['aperture_sum_err'].data[0]]

# 1. HST
SFR_disk_hst, e_SFR_disk_hst = cal_sfr_disk("HST_boundary_1sig_transformed.reg")
print(f"SFR disk (HST) : {SFR_disk_hst:.2f} +/- {e_SFR_disk_hst:.2f} Mo/yr")

# 2. Gemini
SFR_disk_gem, e_SFR_disk_gem = cal_sfr_disk("GMOS_boundary_1sig.reg")
print(f"SFR disk (Gemini) : {SFR_disk_gem:.2f} +/- {e_SFR_disk_gem:.2f} Mo/yr")


# ----- START: [NII]6584 / H alpha flux ratio map ----- #
plt_Data = NII6584_flux_2D / Halpha_flux_2D
plt_Data[plt_Data == 0.] = np.nan
plt_Data[np.isinf(plt_Data) == True] = np.nan

snr_cnd = ((Halpha_snr_2D < 3.0) | (Halpha_snrpix_2D < 5.0) | (NII6584_snr_2D < 3.0))
sig_cnd = (Halpha_sigma_2D < lsig_llim)
rchisq_cnd = (Halpha_rchisq_2D > 50.)
flx_cnd = (Halpha_flux_2D < flx25)
zero_cnd = (snr_cnd | rchisq_cnd)

plt_Data[zero_cnd] = np.nan

fig, ax = plt.subplots(1, 1, figsize=(8,5))
plt.suptitle(r"${\rm [NII]\lambda 6584/H\alpha}$ flux ratio map",
             x=0.5, ha='center', y=0.96, va='top',
             fontsize=20.0)
ax.set_xlim([-3.4, 3.4])
ax.set_ylim([-2.45, 2.45])
ax.set_xticks([-3,-2,-1,0,1,2,3])
ax.set_yticks([-2,-1,0,1,2])
ax.set_xticklabels([r'$-3$',r'$-2$',r'$-1$',0,1,2,3], fontsize=15.0)
ax.set_yticklabels([r'$-2$',r'$-1$',0,1,2], fontsize=15.0)
ax.set_xlabel('arcsec', fontsize=15.0) 
ax.set_ylabel('arcsec', fontsize=15.0)
ax.tick_params(width=1.0, length=5.0)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
im = ax.imshow(plt_Data, cmap='rainbow',
               vmin=0.9*v_low, vmax=1.1*v_high, 
               aspect='equal', extent=[-3.4,3.4,-2.45,2.45])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(im, cax=cax)

ax.contour(X_coord, Y_coord[::-1], sflx, levels=lvs, linewidths=lws, colors=cs, alpha=0.6)
p0, = ax.plot(-100.0, -100.0, '-', linewidth=2.5, color='gray', alpha=0.6,
              label=r"H${\rm \alpha}$ flux contour")
ax.legend(handles=[p0], fontsize=13.0, loc='lower left',
          handlelength=2.5, frameon=True, borderpad=0.8,
          framealpha=0.8, edgecolor='gray')

# The orientations
x0 = -2.75 ; y0 = 1.25
L = 0.6 ; theta0 = gpa*(np.pi/180.0)
ax.arrow(x0-0.025, y0, L*np.sin(theta0), L*np.cos(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.arrow(x0, y0-0.025, -L*np.cos(theta0), L*np.sin(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(-2.95, 2.10, 'E', fontsize=15.0, fontweight='bold', color='blueviolet')
ax.text(-1.90, 1.25, 'N', fontsize=15.0, fontweight='bold', color='blueviolet')

# Scale bar
kpc5 = 5.0 / ang_scale
ax.arrow(2.0, -1.85, kpc5, 0., width=0.07, head_width=0., head_length=0.,
          fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(2.1, -2.2, '5 kpc', fontsize=15.0, fontweight='bold', color='blueviolet')

plt.savefig(dir_fig+'Line_ratio_N2Ha.pdf')
plt.savefig(dir_fig+'Line_ratio_N2Ha.png', dpi=300)
plt.close()
# ----- END: [NII]6584 / H alpha flux ratio map ----- #


# ----- START: Oxygen abundance map ----- #
N2 = np.log10(plt_Data)
logOH = 8.743 + 0.462*N2
plt_Data = logOH

fig, ax = plt.subplots(1, 1, figsize=(8,5))
plt.suptitle(r"${\rm 12+log(O/H)}$ map",
             x=0.5, ha='center', y=0.96, va='top',
             fontsize=20.0)
ax.set_xlim([-3.4, 3.4])
ax.set_ylim([-2.45, 2.45])
ax.set_xticks([-3,-2,-1,0,1,2,3])
ax.set_yticks([-2,-1,0,1,2])
ax.set_xticklabels([r'$-3$',r'$-2$',r'$-1$',0,1,2,3], fontsize=15.0)
ax.set_yticklabels([r'$-2$',r'$-1$',0,1,2], fontsize=15.0)
ax.set_xlabel('arcsec', fontsize=15.0) 
ax.set_ylabel('arcsec', fontsize=15.0)
ax.tick_params(width=1.0, length=5.0)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
im = ax.imshow(plt_Data, cmap='rainbow',
               vmin=8.2, vmax=8.6, 
               aspect='equal', extent=[-3.4,3.4,-2.45,2.45])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(im, cax=cax)

ax.contour(X_coord, Y_coord[::-1], sflx, levels=lvs, linewidths=lws, colors=cs, alpha=0.6)
p0, = ax.plot(-100.0, -100.0, '-', linewidth=2.5, color='gray', alpha=0.6,
              label=r"H${\rm \alpha}$ flux contour")
ax.legend(handles=[p0], fontsize=13.0, loc='lower left',
          handlelength=2.5, frameon=True, borderpad=0.8,
          framealpha=0.8, edgecolor='gray')

# The orientations
x0 = -2.75 ; y0 = 1.25
L = 0.6 ; theta0 = gpa*(np.pi/180.0)
ax.arrow(x0-0.025, y0, L*np.sin(theta0), L*np.cos(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.arrow(x0, y0-0.025, -L*np.cos(theta0), L*np.sin(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(-2.95, 2.10, 'E', fontsize=15.0, fontweight='bold', color='blueviolet')
ax.text(-1.90, 1.25, 'N', fontsize=15.0, fontweight='bold', color='blueviolet')

# Scale bar
kpc5 = 5.0 / ang_scale
ax.arrow(2.0, -1.85, kpc5, 0., width=0.07, head_width=0., head_length=0.,
          fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(2.1, -2.2, '5 kpc', fontsize=15.0, fontweight='bold', color='blueviolet')

plt.savefig(dir_fig+'Metallicity_logOH.pdf')
plt.savefig(dir_fig+'Metallicity_logOH.png', dpi=300)
plt.close()

def cal_fwm_logOH(regfile, data):
    metal = copy.deepcopy(data)
    metal[np.isnan(metal) == True] = 0.

    x0, y0, a, b, pa = read_region(regfile, regtype='ellipse')
    x0, y0 = x0[0]-1, y0[0]-1
    a, b = a[0], b[0]
    pa = pa[0]*np.pi/180.

    ap = EAp(positions=[(x0, y0)], a=a, b=b, theta=pa)

    tbl0 = apphot(data=flx_Data, apertures=ap, mask=None)
    flx_sum_disk = tbl0['aperture_sum'].data[0]
    flx_sum_tail = np.sum(flx_Data) - flx_sum_disk

    tbl1 = apphot(data=flx_Data*metal, apertures=ap, mask=None)
    logOH_fwm_disk = tbl1['aperture_sum'].data[0] / flx_sum_disk
    logOH_fwm_tail = (np.sum(flx_Data*metal) - tbl1['aperture_sum'].data[0]) / flx_sum_tail

    # one_arr = np.ones_like(metal)
    # one_arr[metal == 0.] = 0
    # tbl2 = apphot(data=metal, apertures=ap, mask=None)
    # logOH_mean_disk = tbl2['aperture_sum'].data[0] / (np.pi*a*b)
    # logOH_mean_tail = (np.sum(metal) - tbl2['aperture_sum'].data[0]) / (np.sum(one_arr)-np.pi*a*b)

    return [logOH_fwm_disk, logOH_fwm_tail]

# Total
val = (np.isnan(plt_Data) == False)
fwm_logOH = np.average(plt_Data[val], weights=flx_Data[val])
print(f"Flux-weighted mean of gas metallicity : {fwm_logOH:.3f}")

# 1. HST
logOH_disk_hst, logOH_tail_hst = cal_fwm_logOH("HST_boundary_1sig_transformed.reg", logOH)
print(f"log OH (HST) : {logOH_disk_hst:.3f} + {logOH_tail_hst:.3f}")

# 2. Gemini
logOH_disk_gem, logOH_tail_gem = cal_fwm_logOH("GMOS_boundary_1sig.reg", logOH)
print(f"log OH (Gemini) : {logOH_disk_gem:.3f} + {logOH_tail_gem:.3f}")
# ----- END: Oxygen abundance map ----- #


# ----- START: O3N2 map ----- #
O3N2 = np.log10((OIII5007_flux_2D / Hbeta_flux_2D) * (Halpha_flux_2D / NII6584_flux_2D))
O3N2[np.isinf(O3N2) == True] = np.nan
plt_Data = 10.**O3N2
# plt_Data[np.isinf(plt_Data) == True] = np.nan

fig, ax = plt.subplots(1, 1, figsize=(8,5))
plt.suptitle(r"${\rm ([OIII]\lambda 5007/H\beta)\times(H\alpha/[NII]\lambda5007)}$ map",
             x=0.5, ha='center', y=0.96, va='top',
             fontsize=20.0)
ax.set_xlim([-3.4, 3.4])
ax.set_ylim([-2.45, 2.45])
ax.set_xticks([-3,-2,-1,0,1,2,3])
ax.set_yticks([-2,-1,0,1,2])
ax.set_xticklabels([r'$-3$',r'$-2$',r'$-1$',0,1,2,3], fontsize=15.0)
ax.set_yticklabels([r'$-2$',r'$-1$',0,1,2], fontsize=15.0)
ax.set_xlabel('arcsec', fontsize=15.0) 
ax.set_ylabel('arcsec', fontsize=15.0)
ax.tick_params(width=1.0, length=5.0)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
im = ax.imshow(plt_Data, cmap='rainbow',
               vmin=v_low, vmax=v_high, 
               aspect='equal', extent=[-3.4,3.4,-2.45,2.45])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(im, cax=cax)

ax.contour(X_coord, Y_coord[::-1], sflx, levels=lvs, linewidths=lws, colors=cs, alpha=0.6)
p0, = ax.plot(-100.0, -100.0, '-', linewidth=2.5, color='gray', alpha=0.6,
              label=r"H${\rm \alpha}$ flux contour")
ax.legend(handles=[p0], fontsize=13.0, loc='lower left',
          handlelength=2.5, frameon=True, borderpad=0.8,
          framealpha=0.8, edgecolor='gray')

# The orientations
x0 = -2.75 ; y0 = 1.25
L = 0.6 ; theta0 = gpa*(np.pi/180.0)
ax.arrow(x0-0.025, y0, L*np.sin(theta0), L*np.cos(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.arrow(x0, y0-0.025, -L*np.cos(theta0), L*np.sin(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(-2.95, 2.10, 'E', fontsize=15.0, fontweight='bold', color='blueviolet')
ax.text(-1.90, 1.25, 'N', fontsize=15.0, fontweight='bold', color='blueviolet')

# Scale bar
kpc5 = 5.0 / ang_scale
ax.arrow(2.0, -1.85, kpc5, 0., width=0.07, head_width=0., head_length=0.,
          fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(2.1, -2.2, '5 kpc', fontsize=15.0, fontweight='bold', color='blueviolet')

plt.savefig(dir_fig+'Map_O3N2.pdf')
plt.savefig(dir_fig+'Map_O3N2.png', dpi=300)
plt.close()
# ----- END: O3N2 map ----- #


# ----- START: Oxygen abundance map (2) ----- #
logOH2 = 8.533-0.214*O3N2
plt_Data = logOH2

fig, ax = plt.subplots(1, 1, figsize=(8,5))
plt.suptitle(r"${\rm 12+log(O/H)}$ map",
             x=0.5, ha='center', y=0.96, va='top',
             fontsize=20.0)
ax.set_xlim([-3.4, 3.4])
ax.set_ylim([-2.45, 2.45])
ax.set_xticks([-3,-2,-1,0,1,2,3])
ax.set_yticks([-2,-1,0,1,2])
ax.set_xticklabels([r'$-3$',r'$-2$',r'$-1$',0,1,2,3], fontsize=15.0)
ax.set_yticklabels([r'$-2$',r'$-1$',0,1,2], fontsize=15.0)
ax.set_xlabel('arcsec', fontsize=15.0) 
ax.set_ylabel('arcsec', fontsize=15.0)
ax.tick_params(width=1.0, length=5.0)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
im = ax.imshow(plt_Data, cmap='rainbow',
               vmin=8.2, vmax=8.6, 
               aspect='equal', extent=[-3.4,3.4,-2.45,2.45])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(im, cax=cax)

ax.contour(X_coord, Y_coord[::-1], sflx, levels=lvs, linewidths=lws, colors=cs, alpha=0.6)
p0, = ax.plot(-100.0, -100.0, '-', linewidth=2.5, color='gray', alpha=0.6,
              label=r"H${\rm \alpha}$ flux contour")
ax.legend(handles=[p0], fontsize=13.0, loc='lower left',
          handlelength=2.5, frameon=True, borderpad=0.8,
          framealpha=0.8, edgecolor='gray')

# The orientations
x0 = -2.75 ; y0 = 1.25
L = 0.6 ; theta0 = gpa*(np.pi/180.0)
ax.arrow(x0-0.025, y0, L*np.sin(theta0), L*np.cos(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.arrow(x0, y0-0.025, -L*np.cos(theta0), L*np.sin(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(-2.95, 2.10, 'E', fontsize=15.0, fontweight='bold', color='blueviolet')
ax.text(-1.90, 1.25, 'N', fontsize=15.0, fontweight='bold', color='blueviolet')

# Scale bar
kpc5 = 5.0 / ang_scale
ax.arrow(2.0, -1.85, kpc5, 0., width=0.07, head_width=0., head_length=0.,
          fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(2.1, -2.2, '5 kpc', fontsize=15.0, fontweight='bold', color='blueviolet')

plt.savefig(dir_fig+'Metallicity_logOH2.pdf')
plt.savefig(dir_fig+'Metallicity_logOH2.png', dpi=300)
plt.close()

# Total
val = (np.isnan(plt_Data) == False)
fwm_logOH2 = np.average(plt_Data[val], weights=flx_Data[val])
print(f"Flux-weighted mean of gas metallicity (2) : {fwm_logOH2:.3f}")

# 1. HST
logOH2_disk_hst, logOH2_tail_hst = cal_fwm_logOH("HST_boundary_1sig_transformed.reg", logOH2)
print(f"log OH (HST) : {logOH2_disk_hst:.3f} + {logOH2_tail_hst:.3f}")

# 2. Gemini
logOH2_disk_gem, logOH2_tail_gem = cal_fwm_logOH("GMOS_boundary_1sig.reg", logOH2)
print(f"log OH (Gemini) : {logOH2_disk_gem:.3f} + {logOH2_tail_gem:.3f}")
# ----- END: Oxygen abundance map (2) ----- #


# ----- START: [OIII]5007 / H beta flux ratio map ----- #
plt_Data = OIII5007_flux_2D / Hbeta_flux_2D
plt_Data[plt_Data == 0.] = np.nan
plt_Data[np.isinf(plt_Data) == True] = np.nan

snr_cnd = ((Halpha_snr_2D < 3.0) | (Halpha_snrpix_2D < 5.0) | (OIII5007_snr_2D < 3.0))
sig_cnd = (Halpha_sigma_2D < lsig_llim)
rchisq_cnd = (Halpha_rchisq_2D > 50.)
flx_cnd = (Halpha_flux_2D < flx25)
zero_cnd = (snr_cnd | rchisq_cnd)

plt_Data[zero_cnd] = np.nan

fig, ax = plt.subplots(1, 1, figsize=(8,5))
plt.suptitle(r"${\rm [OIII]\lambda 5007/H\beta}$ flux ratio map",
             x=0.5, ha='center', y=0.96, va='top',
             fontsize=20.0)
ax.set_xlim([-3.4, 3.4])
ax.set_ylim([-2.45, 2.45])
ax.set_xticks([-3,-2,-1,0,1,2,3])
ax.set_yticks([-2,-1,0,1,2])
ax.set_xticklabels([r'$-3$',r'$-2$',r'$-1$',0,1,2,3], fontsize=15.0)
ax.set_yticklabels([r'$-2$',r'$-1$',0,1,2], fontsize=15.0)
ax.set_xlabel('arcsec', fontsize=15.0) 
ax.set_ylabel('arcsec', fontsize=15.0)
ax.tick_params(width=1.0, length=5.0)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
im = ax.imshow(plt_Data, cmap='rainbow',
               vmin=0.9*v_low, vmax=1.1*v_high, 
               aspect='equal', extent=[-3.4,3.4,-2.45,2.45])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(im, cax=cax)

ax.contour(X_coord, Y_coord[::-1], sflx, levels=lvs, linewidths=lws, colors=cs, alpha=0.6)
p0, = ax.plot(-100.0, -100.0, '-', linewidth=2.5, color='gray', alpha=0.6,
              label=r"H${\rm \alpha}$ flux contour")
ax.legend(handles=[p0], fontsize=13.0, loc='lower left',
          handlelength=2.5, frameon=True, borderpad=0.8,
          framealpha=0.8, edgecolor='gray')

# The orientations
x0 = -2.75 ; y0 = 1.25
L = 0.6 ; theta0 = gpa*(np.pi/180.0)
ax.arrow(x0-0.025, y0, L*np.sin(theta0), L*np.cos(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.arrow(x0, y0-0.025, -L*np.cos(theta0), L*np.sin(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(-2.95, 2.10, 'E', fontsize=15.0, fontweight='bold', color='blueviolet')
ax.text(-1.90, 1.25, 'N', fontsize=15.0, fontweight='bold', color='blueviolet')

# Scale bar
kpc5 = 5.0 / ang_scale
ax.arrow(2.0, -1.85, kpc5, 0., width=0.07, head_width=0., head_length=0.,
          fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(2.1, -2.2, '5 kpc', fontsize=15.0, fontweight='bold', color='blueviolet')

plt.savefig(dir_fig+'Line_ratio_O3Hb.pdf')
plt.savefig(dir_fig+'Line_ratio_O3Hb.png', dpi=300)
plt.close()
# ----- END: [OIII]5007 / H beta flux ratio map ----- #


# ----- START: [SII] / H alpha flux ratio map ----- #
plt_Data = (SII6717_flux_2D + SII6731_flux_2D) / Halpha_flux_2D
plt_Data[plt_Data == 0.] = np.nan
plt_Data[np.isinf(plt_Data) == True] = np.nan

snr_cnd = ((Halpha_snr_2D < 3.0) | (Halpha_snrpix_2D < 5.0) | (SII6717_snr_2D < 3.0) | (SII6731_snr_2D < 3.0))
sig_cnd = (Halpha_sigma_2D < lsig_llim)
rchisq_cnd = (Halpha_rchisq_2D > 50.)
flx_cnd = (Halpha_flux_2D < flx25)
zero_cnd = (snr_cnd | rchisq_cnd)

plt_Data[zero_cnd] = np.nan

fig, ax = plt.subplots(1, 1, figsize=(8,5))
plt.suptitle(r"${\rm [SII]\lambda\lambda 6717,6731/H\alpha}$ flux ratio map",
             x=0.5, ha='center', y=0.96, va='top',
             fontsize=20.0)
ax.set_xlim([-3.4, 3.4])
ax.set_ylim([-2.45, 2.45])
ax.set_xticks([-3,-2,-1,0,1,2,3])
ax.set_yticks([-2,-1,0,1,2])
ax.set_xticklabels([r'$-3$',r'$-2$',r'$-1$',0,1,2,3], fontsize=15.0)
ax.set_yticklabels([r'$-2$',r'$-1$',0,1,2], fontsize=15.0)
ax.set_xlabel('arcsec', fontsize=15.0) 
ax.set_ylabel('arcsec', fontsize=15.0)
ax.tick_params(width=1.0, length=5.0)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
im = ax.imshow(plt_Data, cmap='rainbow',
               vmin=0.9*v_low, vmax=1.1*v_high, 
               aspect='equal', extent=[-3.4,3.4,-2.45,2.45])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(im, cax=cax)

ax.contour(X_coord, Y_coord[::-1], sflx, levels=lvs, linewidths=lws, colors=cs, alpha=0.6)
p0, = ax.plot(-100.0, -100.0, '-', linewidth=2.5, color='gray', alpha=0.6,
              label=r"H${\rm \alpha}$ flux contour")
ax.legend(handles=[p0], fontsize=13.0, loc='lower left',
          handlelength=2.5, frameon=True, borderpad=0.8,
          framealpha=0.8, edgecolor='gray')

# The orientations
x0 = -2.75 ; y0 = 1.25
L = 0.6 ; theta0 = gpa*(np.pi/180.0)
ax.arrow(x0-0.025, y0, L*np.sin(theta0), L*np.cos(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.arrow(x0, y0-0.025, -L*np.cos(theta0), L*np.sin(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(-2.95, 2.10, 'E', fontsize=15.0, fontweight='bold', color='blueviolet')
ax.text(-1.90, 1.25, 'N', fontsize=15.0, fontweight='bold', color='blueviolet')

# Scale bar
kpc5 = 5.0 / ang_scale
ax.arrow(2.0, -1.85, kpc5, 0., width=0.07, head_width=0., head_length=0.,
          fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(2.1, -2.2, '5 kpc', fontsize=15.0, fontweight='bold', color='blueviolet')

plt.savefig(dir_fig+'Line_ratio_S2Ha.pdf')
plt.savefig(dir_fig+'Line_ratio_S2Ha.png', dpi=300)
plt.close()
# ----- END: [SII]6731 / H alpha flux ratio map ----- #


# ----- START: [SII]6717 / [SII]6731 flux ratio map ----- #
plt_Data = SII6717_flux_2D / SII6731_flux_2D
plt_Data[plt_Data == 0.] = np.nan
plt_Data[np.isinf(plt_Data) == True] = np.nan

snr_cnd = ((SII6717_snr_2D < 3.0) | (SII6731_snr_2D < 3.0))
sig_cnd = (Halpha_sigma_2D < lsig_llim)
rchisq_cnd = (Halpha_rchisq_2D > 50.)
flx_cnd = (Halpha_flux_2D < flx25)
zero_cnd = (snr_cnd | rchisq_cnd)

plt_Data[zero_cnd] = np.nan

fig, ax = plt.subplots(1, 1, figsize=(8,5))
plt.suptitle(r"${\rm [SII]\lambda 6717/[SII]\lambda 6731}$ flux ratio map",
             x=0.5, ha='center', y=0.96, va='top',
             fontsize=20.0)
ax.set_xlim([-3.4, 3.4])
ax.set_ylim([-2.45, 2.45])
ax.set_xticks([-3,-2,-1,0,1,2,3])
ax.set_yticks([-2,-1,0,1,2])
ax.set_xticklabels([r'$-3$',r'$-2$',r'$-1$',0,1,2,3], fontsize=15.0)
ax.set_yticklabels([r'$-2$',r'$-1$',0,1,2], fontsize=15.0)
ax.set_xlabel('arcsec', fontsize=15.0) 
ax.set_ylabel('arcsec', fontsize=15.0)
ax.tick_params(width=1.0, length=5.0)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
im = ax.imshow(plt_Data, cmap='rainbow',
               vmin=0.9*v_low, vmax=1.1*v_high, 
               aspect='equal', extent=[-3.4,3.4,-2.45,2.45])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(im, cax=cax)

ax.contour(X_coord, Y_coord[::-1], sflx, levels=lvs, linewidths=lws, colors=cs, alpha=0.6)
p0, = ax.plot(-100.0, -100.0, '-', linewidth=2.5, color='gray', alpha=0.6,
              label=r"H${\rm \alpha}$ flux contour")
ax.legend(handles=[p0], fontsize=13.0, loc='lower left',
          handlelength=2.5, frameon=True, borderpad=0.8,
          framealpha=0.8, edgecolor='gray')

# The orientations
x0 = -2.75 ; y0 = 1.25
L = 0.6 ; theta0 = gpa*(np.pi/180.0)
ax.arrow(x0-0.025, y0, L*np.sin(theta0), L*np.cos(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.arrow(x0, y0-0.025, -L*np.cos(theta0), L*np.sin(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(-2.95, 2.10, 'E', fontsize=15.0, fontweight='bold', color='blueviolet')
ax.text(-1.90, 1.25, 'N', fontsize=15.0, fontweight='bold', color='blueviolet')

# Scale bar
kpc5 = 5.0 / ang_scale
ax.arrow(2.0, -1.85, kpc5, 0., width=0.07, head_width=0., head_length=0.,
          fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(2.1, -2.2, '5 kpc', fontsize=15.0, fontweight='bold', color='blueviolet')

plt.savefig(dir_fig+'Line_ratio_S2S2.pdf')
plt.savefig(dir_fig+'Line_ratio_S2S2.png', dpi=300)
plt.close()
# ----- END: [SII]6717 / [SII]6731 flux ratio map ----- #


# ----- START: [NII]6584 / [NII]6548 flux ratio map ----- #
plt_Data = NII6584_flux_2D / NII6548_flux_2D
plt_Data[plt_Data == 0.] = np.nan
plt_Data[np.isinf(plt_Data) == True] = np.nan

snr_cnd = ((NII6548_snr_2D < 3.0) | (NII6584_snr_2D < 3.0))
sig_cnd = (Halpha_sigma_2D < lsig_llim)
rchisq_cnd = (Halpha_rchisq_2D > 50.)
flx_cnd = (Halpha_flux_2D < flx25)
zero_cnd = (snr_cnd | rchisq_cnd)

plt_Data[zero_cnd] = np.nan

fig, ax = plt.subplots(1, 1, figsize=(8,5))
plt.suptitle(r"${\rm [NII]\lambda 6584/[NII]\lambda 6548}$ flux ratio map",
             x=0.5, ha='center', y=0.96, va='top',
             fontsize=20.0)
ax.set_xlim([-3.4, 3.4])
ax.set_ylim([-2.45, 2.45])
ax.set_xticks([-3,-2,-1,0,1,2,3])
ax.set_yticks([-2,-1,0,1,2])
ax.set_xticklabels([r'$-3$',r'$-2$',r'$-1$',0,1,2,3], fontsize=15.0)
ax.set_yticklabels([r'$-2$',r'$-1$',0,1,2], fontsize=15.0)
ax.set_xlabel('arcsec', fontsize=15.0) 
ax.set_ylabel('arcsec', fontsize=15.0)
ax.tick_params(width=1.0, length=5.0)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
im = ax.imshow(plt_Data, cmap='rainbow',
               vmin=0.9*v_low, vmax=1.1*v_high, 
               aspect='equal', extent=[-3.4,3.4,-2.45,2.45])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(im, cax=cax)

ax.contour(X_coord, Y_coord[::-1], sflx, levels=lvs, linewidths=lws, colors=cs, alpha=0.6)
p0, = ax.plot(-100.0, -100.0, '-', linewidth=2.5, color='gray', alpha=0.6,
              label=r"H${\rm \alpha}$ flux contour")
ax.legend(handles=[p0], fontsize=13.0, loc='lower left',
          handlelength=2.5, frameon=True, borderpad=0.8,
          framealpha=0.8, edgecolor='gray')

# The orientations
x0 = -2.75 ; y0 = 1.25
L = 0.6 ; theta0 = gpa*(np.pi/180.0)
ax.arrow(x0-0.025, y0, L*np.sin(theta0), L*np.cos(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.arrow(x0, y0-0.025, -L*np.cos(theta0), L*np.sin(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(-2.95, 2.10, 'E', fontsize=15.0, fontweight='bold', color='blueviolet')
ax.text(-1.90, 1.25, 'N', fontsize=15.0, fontweight='bold', color='blueviolet')

# Scale bar
kpc5 = 5.0 / ang_scale
ax.arrow(2.0, -1.85, kpc5, 0., width=0.07, head_width=0., head_length=0.,
          fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(2.1, -2.2, '5 kpc', fontsize=15.0, fontweight='bold', color='blueviolet')

plt.savefig(dir_fig+'Line_ratio_N2N2.pdf')
plt.savefig(dir_fig+'Line_ratio_N2N2.png', dpi=300)
plt.close()
# ----- END: [NII]6584 / [NII]6548 flux ratio map ----- #


# ----- START: BPT diagram ----- #
x_Dat = np.log10(NII6584_flux_2D/Halpha_flux_2D)
y_Dat = np.log10(OIII5007_flux_2D/Hbeta_flux_2D)

snr_cnd = ((NII6584_snr_2D < 3.0) | (Halpha_snr_2D < 3.0) | \
           (OIII5007_snr_2D < 3.0) | (Hbeta_snr_2D < 3.0))
sig_cnd = (Halpha_sigma_2D < lsig_llim)
rchisq_cnd = (Halpha_rchisq_2D > 50.)
flx_cnd = (Halpha_flux_2D < flx25)
zero_cnd = (snr_cnd | rchisq_cnd)

x_Dat[zero_cnd] = np.nan
y_Dat[zero_cnd] = np.nan

BPT_SFG = (y_Dat < 1.3+0.61/x_Dat)
BPT_comp = ((y_Dat > 1.3+0.61/x_Dat) & (y_Dat < 1.19+0.61/(x_Dat-0.47)))
BPT_AGN = ((y_Dat > 1.19+0.61/(x_Dat-0.47)) & \
           (y_Dat-1.4899 > ((0.2949- -0.2135)/(1.4899-0.4548))*(x_Dat-0.2949)))
BPT_LINER = ((y_Dat > 1.19+0.61/(x_Dat-0.47)) & \
             (y_Dat-1.4899 < ((0.2949- -0.2135)/(1.4899-0.4548))*(x_Dat-0.2949)))

x_Dat_flat = x_Dat.flatten()
y_Dat_flat = y_Dat.flatten()

fig, ax = plt.subplots(1, 1, figsize=(8,5))
plt.suptitle("BPT diagram",
             x=0.5, ha='center', y=0.96, va='top',
             fontsize=20.0)
ax.set_xlim([-1.4, 0.5])
ax.set_ylim([-1.0, 1.5])
ax.set_xlabel(r'log([NII]$\lambda$6584/H$\alpha$)', fontsize=15.0) 
ax.set_ylabel(r'log([OIII]$\lambda$5007/H$\beta$)', fontsize=15.0)
ax.tick_params(axis='both', labelsize=15.0)
ax.tick_params(width=1.0, length=5.0)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

BPT_class = ['SFG','comp','AGN','LINER']
BPT_color = ['dodgerblue','darkorange','red','green']
BPT_cmap = ['Blues', 'Oranges', 'Reds', 'Greens']

for i in np.arange(len(BPT_class)):
    exec("X = x_Dat_flat[BPT_"+BPT_class[i]+".flatten()]")
    exec("Y = y_Dat_flat[BPT_"+BPT_class[i]+".flatten()]")
    uniq_val, uniq_idx = np.unique(X, return_index=True)
    if (len(uniq_idx) > 0):
        cmap = cm.get_cmap(BPT_cmap[i])
        c_lo, c_hi = np.min(Y[uniq_idx])-0.1, np.max(Y[uniq_idx])+0.1
        for j in np.arange(len(uniq_idx)):
            c_id = 0.1+(0.9-0.1)*(Y[uniq_idx][j]-c_lo)/(c_hi-c_lo)
            ax.plot(X[uniq_idx][j], Y[uniq_idx][j],
                    'o', ms=5.0, color=cmap(c_id)[:-1], alpha=0.8)

# Kauffmann+03
xx = np.linspace(-1.5, -0.01, 1000)
ax.plot(xx, 1.3+0.61/xx, 'k-', linewidth=1.25, alpha=0.75)

# Kewley+01
xx = np.linspace(-1.5, 0.45, 1000)
ax.plot(xx, 1.19+0.61/(xx-0.47), 'k:', linewidth=1.25, alpha=0.75)

# AGN-LINER boundary (digitized from GASP XV. paper)
ax.plot([-0.2135, 0.4548], [0.2949, 1.4899], 'k--', linewidth=1.25, alpha=0.75)

plt.savefig(dir_fig+'BPT_diagram.pdf')
plt.savefig(dir_fig+'BPT_diagram.png', dpi=300)
plt.close()
# ----- END: BPT diagram ----- #


# ----- START: BPT spatial map ----- #
fig, ax = plt.subplots(1, 1, figsize=(8,5))
plt.suptitle("The spatial map of the BPT class",
             x=0.5, ha='center', y=0.96, va='top',
             fontsize=20.0)
ax.set_xlim([-3.4, 3.4])
ax.set_ylim([-2.45, 2.45])
ax.set_xticks([-3,-2,-1,0,1,2,3])
ax.set_yticks([-2,-1,0,1,2])
ax.set_xticklabels([r'$-3$',r'$-2$',r'$-1$',0,1,2,3], fontsize=15.0)
ax.set_yticklabels([r'$-2$',r'$-1$',0,1,2], fontsize=15.0)
ax.set_xlabel('arcsec', fontsize=15.0) 
ax.set_ylabel('arcsec', fontsize=15.0)
ax.tick_params(width=1.0, length=5.0)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

for i in np.arange(len(BPT_class)):
    exec("bpt_data = BPT_"+BPT_class[i])
    if (np.sum(bpt_data == True) > 0):
        cmap = cm.get_cmap(BPT_cmap[i])
        c_lo, c_hi = np.min(y_Dat[bpt_data == True])-0.1, np.max(y_Dat[bpt_data == True])+0.1
        for j in np.arange(np.sum(bpt_data == True)):
            c_id = 0.1+(0.9-0.1)*(y_Dat[bpt_data == True][j]-c_lo)/(c_hi-c_lo)
            rect = plt.Rectangle((-0.1*bpt_data.shape[1]/2 + 0.1*np.argwhere(bpt_data == True)[j,1],
                                  0.1*bpt_data.shape[0]/2 - 0.1*np.argwhere(bpt_data == True)[j,0]),
                                 0.1, 0.1, facecolor=cmap(c_id)[:-1], alpha=0.8)
            ax.add_patch(rect)

ax.contour(X_coord, Y_coord[::-1], sflx, levels=lvs, linewidths=lws, colors=cs, alpha=0.6)
p0, = ax.plot(-100.0, -100.0, '-', linewidth=2.5, color='gray', alpha=0.6,
              label=r"H${\rm \alpha}$ flux contour")
ax.legend(handles=[p0], fontsize=13.0, loc='lower left',
          handlelength=2.5, frameon=True, borderpad=0.8,
          framealpha=0.8, edgecolor='gray')

# The orientations
x0 = -2.75 ; y0 = 1.25
L = 0.6 ; theta0 = gpa*(np.pi/180.0)
ax.arrow(x0-0.025, y0, L*np.sin(theta0), L*np.cos(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.arrow(x0, y0-0.025, -L*np.cos(theta0), L*np.sin(theta0), width=0.06,
         head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(-2.95, 2.10, 'E', fontsize=15.0, fontweight='bold', color='blueviolet')
ax.text(-1.90, 1.25, 'N', fontsize=15.0, fontweight='bold', color='blueviolet')

# Scale bar
kpc5 = 5.0 / ang_scale
ax.arrow(2.0, -1.85, kpc5, 0., width=0.07, head_width=0., head_length=0.,
          fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(2.1, -2.2, '5 kpc', fontsize=15.0, fontweight='bold', color='blueviolet')

plt.savefig(dir_fig+'BPT_map.pdf')
plt.savefig(dir_fig+'BPT_map.png', dpi=300)
plt.close()
# ----- END: BPT spatial map ----- #


# ----- Saving the results ----- #
np.savez('plot_data.npz', sflx=sflx, rvd=Rvd, vdd=Vdd, sfrd=SFRD)

'''
# ----- Applying NebulaBayes ----- #

# These HII-region optical emission-line fluxes have already been dereddened
linelist = ["Hbeta", "OIII5007", "Halpha", "NII6583", "SII6716", "SII6731"]
obs_fluxes = [8.151, 4.634, 1.999, 0.09562, 13.21, 5.116, 1.377, 1.249]
obs_errs = [0.09008, 0.04013, 0.01888, 0.00222, 0.07635, 0.03159, 0.00999,
            0.00923]
# Fluxes/errors will be normalised to the flux of the default norm_line, Hbeta
# The emission line names match those in the grid (see the NebulaBayes/grids
# directory)

# Initialise the NB_Model, which loads and interpolates the model flux grids:
NB_Model_HII = NB_Model("HII", line_list=linelist)

# Set outputs:
kwargs = {"prior_plot": os.path.join(OUT_DIR, "1_HII_prior_plot.pdf"),
          "likelihood_plot": os.path.join(OUT_DIR, "1_HII_likelihood_plot.pdf"),
          "posterior_plot": os.path.join(OUT_DIR, "1_HII_posterior_plot.pdf"),
          "estimate_table": os.path.join(OUT_DIR, "1_HII_param_estimates.csv"),
          "best_model_table": os.path.join(OUT_DIR, "1_HII_best_model.csv"),
          }

# Run parameter estimation once
Result_HII = NB_Model_HII(obs_fluxes, obs_errs, linelist, **kwargs)
# NB_Model_HII may be called repeatedly to do Bayesian parameter estimation
# on different sets of observed fluxes with the same grid.
'''



# Printing the running time
print('\n')
print('--- %s seconds ---' %(time.time()-start_time))
