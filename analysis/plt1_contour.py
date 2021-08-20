#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 09:59:44 2020

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import copy
from astropy.io import fits
from astropy import wcs
import imgscl
from scipy import ndimage
from skimage.registration import phase_cross_correlation
from astropy.cosmology import FlatLambdaCDM


# ----- Loading the original HST image ----- #

# Directories
dir_fig = '/data/jlee/DATA/Gemini/Programs/GN-2019A-Q-215/analysis/diagram/linefits/'
diI = '/data/jlee/DATA/HLA/McPartland+16/MACS1752/JFG2/Phot/'
diG = '/data/jlee/DATA/Gemini/Programs/GN-2019A-Q-215/redux4_700/'

# Reading central RA & Dec
hdr1 = fits.getheader(diG+'cstxeqxbrgN20190611S0257_3D.fits', ext=0)
gra, gdec, gpa = hdr1['RA'], hdr1['DEC'], hdr1['PA']

# Reading FITS images and creating RGB data
#img435 = fits.getdata(diI+'alw435NE_drc.fits')
img606 = fits.getdata(diI+'606_ori.fits', ext=0)
img814 = fits.getdata(diI+'814_ori.fits', ext=0)

cimg = np.zeros((img814.shape[0], img814.shape[1], 3), dtype=float)
cimg[:,:,0] = imgscl.linear(0.5*img814, scale_min=-0.02, scale_max=0.075)   # R
cimg[:,:,1] = imgscl.linear(0.5*(img606+img814), scale_min=-0.02, scale_max=0.15)   # G
cimg[:,:,2] = imgscl.linear(0.5*img606, scale_min=-0.02, scale_max=0.075)   # B

# WCS to XY
h = fits.getheader(diI+'606_ori.fits', ext=0)
w = wcs.WCS(h)
px, py = w.wcs_world2pix(gra, gdec, 1)


# Color map
rth = 150.0
img = cimg[int(np.round(py)-1-rth):int(np.round(py)-1+rth),
           int(np.round(px)-1-rth):int(np.round(px)-1+rth),:]


# ----- Rotating the image ----- #
rotated_img = ndimage.rotate(img, gpa)
# img -> rotated_img


# ----- Background values ----- #
m0, m1, m2 = np.nanmean(img[210:290, 20:100], axis=(0,1))
s0, s1, s2 = np.std(img[210:290, 20:100], axis=(0,1))
# print(m0, m1, m2)
# print(s0, s1, s2)


# ----- Putting pixel values to no-signal pixels ----- #
rimg = copy.deepcopy(rotated_img)
no_sign = (rimg == [0., 0., 0.])
no_sign_idx = np.argwhere(no_sign)

for i in np.arange(no_sign_idx.shape[0]):
    rimg[no_sign_idx[i,0], no_sign_idx[i,1], :] = np.random.normal([m0, m1, m2], [s0, s1, s2])
# img -> rotated_img -> rimg


# ----- Reading H alpha data ----- #
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
redshift = 0.3527
dist_lum = cosmo.luminosity_distance(redshift).value*1.0e+6    # pc
c = 2.99792e+5    # km/s
ang_scale = cosmo.kpc_proper_per_arcmin(redshift).value / 60.    # kpc/arcsec
pixel_scale = 0.1    # arcsec/pix

wav_Ha = 6564.61  # Angstrom
par, e_par = np.loadtxt('relation_wav_R.txt').T
Rmax = par[0] + par[1]*wav_Ha*(1+redshift)
e_Rmax = np.sqrt(e_par[0]**2.0 + (e_par[1]*wav_Ha*(1+redshift))**2.0)
lsig0 = wav_Ha / (2.0*np.sqrt(2.0*np.log(2.0))*Rmax)
e_lsig0 = wav_Ha*e_Rmax / (2.0*np.sqrt(2.0*np.log(2.0))*Rmax*Rmax)
lsig_llim = lsig0 - 1.0*e_lsig0
vsig0 = c / (2.0*np.sqrt(2.0*np.log(2.0))*Rmax)

dir_Ha = '/data/jlee/DATA/Gemini/Programs/GN-2019A-Q-215/analysis/lines2/Halpha/'
Ha_flx = fits.getdata(dir_Ha+'flux_2D.fits', ext=0)
Ha_snr = fits.getdata(dir_Ha+'snr_2D.fits', ext=0)
Ha_sig = fits.getdata(dir_Ha+'sigma_2D.fits', ext=0)
Ha_rchisq = fits.getdata(dir_Ha+'rchisq_2D.fits', ext=0)
Ha_snrpix = fits.getdata(dir_Ha+'snrpix_2D.fits', ext=0)

snr_cnd = ((Ha_snr < 3.0) | (Ha_snrpix < 5.0))
# sig_cnd = (Ha_sig < lsig_llim)
rchi25, rchi50, rchi75 = np.percentile(Ha_rchisq[Ha_rchisq > 0.], [25.0, 50.0, 75.0])
rchisq_cnd = (Ha_rchisq > 50.)
flx25, flx50, flx75 = np.percentile(Ha_flx[Ha_flx > 0.], [25.0, 50.0, 75.0])
flx_cnd = (Ha_flx < flx50)
zero_cnd = (snr_cnd | rchisq_cnd)# | flx_cnd)
# zero_cnd = (snr_cnd | flx_cnd)
Ha_flx[zero_cnd] = 0.
sflx = ndimage.gaussian_filter(Ha_flx, sigma=(1.25, 1.25), order=0)


# ----- Making cutout image for 2x2 binning ----- #
hstimg = rimg[int(0.5*(rimg.shape[0]-1)-Ha_flx.shape[0]):int(0.5*(rimg.shape[0]-1)+Ha_flx.shape[0]),
              int(0.5*(rimg.shape[1]-1)-Ha_flx.shape[1]):int(0.5*(rimg.shape[1]-1)+Ha_flx.shape[1]), :]
# img -> rotated_img -> rimg -> hstimg


# ----- 2x2 binning for matching the pixel scale ----- #
bin_nx, bin_ny = 2, 2
hstimg_binned = np.zeros((hstimg.shape[0]//bin_ny,
                          hstimg.shape[1]//bin_nx,
                          hstimg.shape[2]))
for y in np.arange(int(hstimg.shape[0]/bin_ny)):
    for x in np.arange(int(hstimg.shape[1]/bin_nx)):
        hstimg_binned[y, x, :] = np.mean(hstimg[bin_ny*y:bin_ny*(y+1), bin_nx*x:bin_nx*(x+1), :], axis=(0,1))
# img -> rotated_img -> rimg -> hstimg -> hstimg_binned


# ----- Making the 2D image ----- #
hstimg_2D = np.sum(hstimg_binned, axis=2)


# ----- Running the 2D correlation ----- #
shifted, error, diffphase = phase_cross_correlation(np.flip(Ha_flx, axis=0), hstimg_2D, upsample_factor=100)
print(shifted)
corr_hstimg = ndimage.shift(hstimg, shift=(2*shifted[0], 2*shifted[1], 0), mode='nearest')
Y_coord = pixel_scale*(np.arange(hstimg_2D.shape[0], step=1)-hstimg_2D.shape[0]/2)
X_coord = pixel_scale*(np.arange(hstimg_2D.shape[1], step=1)-hstimg_2D.shape[1]/2)
# img -> rotated_img -> rimg -> hstimg -> hstimg_binned -> corr_hstimg


# ----- Plotting figure ----- #
def plot_contour(hst_Data, sflux_Data, out,
                 x0=-2.75, y0=1.25, L=0.6, theta0=gpa*(np.pi/180.0),
                 xN=-1.90, yN=1.25, xE=-2.95, yE=2.10):

    fig, ax = plt.subplots(1, 1, figsize=(8,5))

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

    # HST image
    im = ax.imshow(hst_Data, extent=[-3.4,3.4,-2.45,2.45],
                   origin='lower', aspect='equal')

    sig = np.std(sflux_Data)
    lvs = [0.5*sig, 1.*sig, 2.*sig, 3.*sig, 5.*sig, 7.*sig]
    lws = tuple(np.repeat(3.75, len(lvs)-1))
    cs = tuple(['magenta']*(len(lvs)-1))

    # GMOS/IFU contour
    ax.contour(X_coord, Y_coord[::-1], sflx, levels=lvs, linewidths=lws, colors=cs, alpha=0.7)
    p0, = ax.plot(-100.0, -100.0, '-', linewidth=lws[0], color=cs[0], alpha=0.7,
                  label=r"H${\rm \alpha}$ flux contour")

    # The orientations
    ax.arrow(x0-0.025, y0, L*np.sin(theta0), L*np.cos(theta0), width=0.06,
             head_width=0.18, head_length=0.18, fc='yellow', ec='yellow', alpha=0.9)
    ax.arrow(x0, y0-0.025, -L*np.cos(theta0), L*np.sin(theta0), width=0.06,
             head_width=0.18, head_length=0.18, fc='yellow', ec='yellow', alpha=0.9)
    ax.text(xE, yE, 'E', fontsize=17.5, fontweight='bold', color='yellow')
    ax.text(xN, yN, 'N', fontsize=17.5, fontweight='bold', color='yellow')

    # Scale bar
    kpc5 = 5.0 / ang_scale
    ax.arrow(2.0, -1.85, kpc5, 0., width=0.07, head_width=0., head_length=0.,
              fc='yellow', ec='yellow', alpha=0.9)
    ax.text(2.1, -2.2, '5 kpc', fontsize=17.5, fontweight='bold', color='yellow')

    # Legend
    ax.legend(handles=[p0], fontsize=13.0, loc='lower left',
              handlelength=2.5, frameon=True, borderpad=0.8,
              framealpha=0.8, edgecolor='gray')

    plt.tight_layout()

    plt.savefig(out+'.pdf', dpi=300)
    plt.savefig(out+'.png', dpi=300)
    plt.close()

plot_contour(corr_hstimg, sflx, dir_fig+"Halpha_contour",
             x0=-2.75, y0=1.25, L=0.6, theta0=gpa*(np.pi/180.0),
             xN=-1.90, yN=1.25, xE=-2.95, yE=2.10)


# ----- Saving the results ----- #
np.savez('contour_data.npz', img=corr_hstimg, sflx=sflx, x=X_coord, y=Y_coord)


# Printing the running time
print('\n')
print('--- %s seconds ---' %(time.time()-start_time))