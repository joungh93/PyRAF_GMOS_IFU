#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 16:20:45 2020

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import glob, os, copy
from matplotlib import pyplot as plt
import init_cfg as ic
from astropy.io import fits
import vorbin
from vorbin.voronoi_2d_binning import voronoi_2d_binning
from scipy.stats import sigmaclip
from scipy.optimize import curve_fit
from astropy.convolution import convolve
from astropy.convolution import Gaussian1DKernel
from tqdm import trange
from astropy.cosmology import FlatLambdaCDM
from specutils import Spectrum1D
from specutils.fitting import fit_continuum
from astropy import units as u
from astropy.modeling.polynomial import Chebyshev1D
import warnings
warnings.filterwarnings("ignore")


# ----- Basic parameters ----- #
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
redshift = ic.redshift
dist_lum = cosmo.luminosity_distance(redshift).value*1.0e+6    # pc

dir_vbin = 'vorbin/'
dir_vbin0 = 'vorbin_Test/'
os.system('rm -rfv '+dir_vbin)
os.system('mkdir '+dir_vbin)


# ----- Reading the cube ----- #
fin_cb = 'bfcube_3D.fits'

hd0 = fits.getheader(fin_cb, ext=0)
d_sci, h_sci = fits.getdata(fin_cb, ext=1, header=True)
d_var, h_var = fits.getdata(fin_cb, ext=2, header=True)

wav = np.linspace(start=h_sci['CRVAL3']+(1-h_sci['CRPIX3'])*h_sci['CD3_3'],
                  stop=h_sci['CRVAL3']+(h_sci['NAXIS3']-h_sci['CRPIX3'])*h_sci['CD3_3'],
                  num=h_sci['NAXIS3'], endpoint=True)

wav_rest = wav/(1.0+redshift)
# d_sci *= (1.0+redshift)
# d_var *= (1.0+redshift)**2.0

d_sci_sb = copy.deepcopy(d_sci)
d_var_sb = copy.deepcopy(d_var)

d_sci_sb[:, 0, :] = 0.
d_sci_sb[:, :, 0] = 0.
d_sci_sb[:, -1, :] = 0.
d_sci_sb[:, :, -1] = 0.


# ----- Voronoi binning ----- #
# Wavelength range: H alpha + [NII] lines
# d_snr = np.maximum(0, d_sci) / np.sqrt(d_var)
wav_rest_range = [[6545.0, 6590.0]]#, [4855.0, 4870.0]]
wav_range = (1+redshift)*np.array(wav_rest_range)
spx_range = []
for w in np.arange(len(wav_range)):
    spx_range.append([np.abs(wav-wav_range[w,0]).argmin(),
                      np.abs(wav-wav_range[w,1]).argmin()])
spx_range = np.array(spx_range)

# Writing & reloading the input text file
x, y, sig, noi = np.loadtxt(dir_vbin0+'vorbin_input2.txt').T


# Voronoi 2D binning
targetSN = 50.0

fig, ax = plt.subplots()
binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning(
        x, y, sig, noi, targetSN, plot=1, quiet=1)
plt.savefig(dir_vbin+'vbin.png', dpi=300)
plt.close()

uniq_val, uniq_ind, uniq_cnt = np.unique(binNum, return_index=True, return_counts=True) 
print(f'Total number of bin: {len(x):d}')
print(f'After Voronoi binning: {len(uniq_val):d}')
print(f'1-pixel bin ratio: {100.0*np.sum(uniq_cnt == 1)/len(uniq_val):.2f}%')


# Spectra binning
ix, iy = x.astype(int), y.astype(int)
nvbin = sn.size
binned_spec = np.zeros([d_sci_sb.shape[0], nvbin], dtype=float)
binned_vari = np.zeros_like(binned_spec)
for i in np.arange(nvbin):
    idxs = (binNum == i)
    for idx, idy in list(zip(ix[idxs], iy[idxs])):
        var_xy = d_var_sb[:, idy, idx]
        var_xy[((var_xy <= 0.) | (np.isnan(var_xy)))] = np.nanmedian(var_xy)
        binned_spec[:, i] += d_sci_sb[:, idy, idx]
        binned_vari[:, i] += var_xy
bin_order = sn.argsort()[::-1]


# Wavelength masking range
wav_frame_mode = 'obs'    # 'obs' / 'res'
wav_msk2 = np.array([[4990, 5015],  # [OII]
                    [5170, 5200],   # [NeIII]
                    [5430, 5560],   # H delta
                    [5750, 5875],   # H gamma
                    [6000, 6320],   # noisy? region
                    [6450, 6580],   # H beta
                    [6630, 6780],    # [OIII]4959/5007
                    [6900, 7300],   # noisy? region
                    [8700, 8900],   # [NII] + H alpha
                    [9000, 9075]    # [SII]
                    ])
if (wav_frame_mode == 'res'):
    wav_msk2 = wav_msk2*(1.0+redshift)
else:
    pass

binned_cont = np.zeros_like(binned_spec)
for i in trange(nvbin):
    fig, ax = plt.subplots()
    ax.set_xlabel(r'Angstrom [${\rm \AA}$] (Observer-frame)')
    ax.set_ylabel(r'Flux [$10^{-15}~{\rm erg~cm^{-2}~s^{-1}~\AA^{-1}}$]')
    p1, = ax.plot(wav, binned_spec[:, bin_order[i]]/nPixels[bin_order[i]],
                  '-', color='dodgerblue', linewidth=2.0, alpha=0.5, zorder=1,
                  label='Binned spectrum')

    # Continuum fitting: polymial relation
    region = [(wav.min()*u.AA, wav_msk2[0,0]*u.AA)]
    for j in np.arange(wav_msk2.shape[0]-1):
        region.append((wav_msk2[j,1]*u.AA, wav_msk2[j+1,0]*u.AA))
    region.append((wav_msk2[wav_msk2.shape[0]-1, 1]*u.AA, wav.max()*u.AA))
    # print(region)

    lam = wav * u.AA
    flx = binned_spec[:, bin_order[i]]/nPixels[bin_order[i]] * u.Unit("1.0e-15 erg cm-2 s-1 AA-1")
    vspec = Spectrum1D(spectral_axis = lam, flux = flx)
    cfunc = fit_continuum(vspec, model=Chebyshev1D(10), window=region)
    ccont = cfunc(vspec.spectral_axis)
    p4, = ax.plot(wav, ccont.value, '-', color='red', linewidth=1.5, alpha=0.7, zorder=5,
                  label='Chebyshev continuum')

    plt.legend(handles=[p1, p4], fontsize=9.0, loc='upper left',
               handlelength=2.5, frameon=True, borderpad=0.8)
    plt.tight_layout()
    plt.savefig(dir_vbin+f'vbin_{bin_order[i]:03d}.png', dpi=300)
    plt.close()

    binned_cont[:, bin_order[i]] = ccont.value * nPixels[bin_order[i]]


# Saving the results
np.savetxt(dir_vbin+'vorbin_output.txt', np.column_stack([x, y, binNum, sn[binNum], nPixels[binNum], scale[binNum]]),
           fmt='%4i  %4i  %4i  %6.2f  %4i  %6.2f')

data_vbin = np.zeros((d_sci.shape[1], d_sci.shape[2]))
data_vbin[:, :] = -99
for i in np.arange(len(binNum)):
    data_vbin[iy[i], ix[i]] = binNum[i]
fits.writeto(dir_vbin+'vbin.fits', data_vbin, overwrite=True)

np.savez(dir_vbin+'vorbin_array.npz', wav=wav, sci=binned_spec, var=binned_vari, cont=binned_cont)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
