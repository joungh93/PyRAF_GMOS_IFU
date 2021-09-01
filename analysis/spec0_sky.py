#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 12:09:25 2020

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import glob, os
import pandas as pd
import init_cfg as ic
from astropy.io import fits
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
from scipy import interpolate
import copy
from scipy.special import erf
from astropy.stats import sigma_clip
import tqdm
from scipy.odr import *


# ----- Basic parameters & function ----- #
dir_fig = ic.cpath+"diagram/"
c = 2.99792e+5
def gauss_cdf_scale(x, mu, sigma, flux_scale):
    dx = x[1] - x[0]
    v1 = erf((x-mu+0.5*dx)/(np.sqrt(2.0)*sigma))
    v2 = erf((x-mu-0.5*dx)/(np.sqrt(2.0)*sigma))
    return flux_scale*(v1-v2)/(2.0*dx)


# ----- Making a list of sky frames ----- #
sky_list = []
for i in np.arange(len(ic.cube_list)):
    di = ic.cube_list[i].split('cstxeqxbrg')[0]
    cb = ic.cube_name[i]
    sky_list.append(di+'stxeqxbrg'+cb+'.fits')


# ----- Wavelength alignment ------ #
wav_start, wav_end = ic.wav_range[0]+10., ic.wav_range[1]-10.
wav_new = np.linspace(start=wav_start, stop=wav_end,
                      num=1+int((wav_end-wav_start)/ic.wav_intv), endpoint=True)

d_sky = np.zeros((len(sky_list), 1+int((wav_end-wav_start)/ic.wav_intv)))
for i in np.arange(len(sky_list)):
    hdr_sci = fits.getheader(sky_list[i], ext=2)
    dat_sky, hdr_sky = fits.getdata(sky_list[i], ext=5, header=True)
    
    wav = np.linspace(start=hdr_sci['CRVAL1']+(1-hdr_sci['CRPIX1'])*hdr_sci['CD1_1'],
                      stop=hdr_sci['CRVAL1']+(hdr_sci['NAXIS1']-hdr_sci['CRPIX1'])*hdr_sci['CD1_1'],
                      num=hdr_sci['NAXIS1'], endpoint=True) 

    nw = int(round((ic.wav_range[1]-ic.wav_range[0])/hdr_sci['CD1_1']))
    spx_start = np.abs(wav-ic.wav_range[0]).argmin()
    spx_end = spx_start + nw
    wav_input = wav[spx_start:spx_end]

    func_sky = interpolate.interp1d(wav_input, dat_sky[spx_start:spx_end], kind='linear')
    dat_sky_new = func_sky(wav_new)

    d_sky[i, :] = dat_sky_new

if (ic.combine_mode == 'median'):
    cd_sky = np.nanmedian(d_sky, axis=0)
if (ic.combine_mode == 'clippedmean'):
    clipped_ds = sigma_clip(d_sky, sigma=3.0, axis=0).data
    cd_sky = np.nanmean(clipped_ds, axis=0)

d_sky = cd_sky

col1 = fits.Column(name='wavelength', format='E', array=wav_new)
col2 = fits.Column(name='count', format='E', array=d_sky)
h = fits.BinTableHDU.from_columns([col1, col2])
h.writeto('sky.fits', overwrite=True)


### Figure 1 : sky spectra
fig = plt.figure(1, figsize=(8,4))
ax = plt.subplot(1,1,1)
ax.set_position([0.10, 0.17, 0.85, 0.77])
ax.set_xlim([4500.0, 10000.0])
ax.set_ylim([-50.0, 1000.0])
ax.set_xticks(np.linspace(5000, 10000, 6))
ax.set_yticks(np.linspace(0, 1000, 6))
ax.set_xticklabels(['5000','6000','7000','8000','9000','10000'], fontsize=13)
ax.set_yticklabels(['0','200','400','600','800','1000'], fontsize=13)
ax.set_xlabel(r'Wavelength $[\AA]$', fontsize=13)
ax.set_ylabel('Count', fontsize=13)
ax.tick_params(width=1.5, length=8.0)
plt.minorticks_on()
ax.tick_params(width=1.5,length=5.0,which='minor')
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)

ax.plot(wav_new, d_sky, '-', color='tab:blue', linewidth=1.5, alpha=0.8)

plt.savefig(dir_fig+'sky_spectra.png', dpi=300)
plt.close()


# ----- Sky line fitting ----- #

# Fitting 1 : 5570-5585AA
# Fitting 2 : OH 7-2 Q1 (1.5) - 6863.955AA
# Fitting 3 : 8390-8405AA
# Fitting 4 : 8875-8895AA
# Fitting 5 : 8895-8910AA
# Fitting 6 : 8910-8930AA

wav_fit = [[5570., 5585.],  # a skyline between 5570-5585AA
           [6855., 6870.],  # OH 7-2 Q1 (1.5) - 6863.955AA
           [8390., 8405.],  # a skyline between 8390-8405AA
           [8880., 8895.],  # a skyline between 8875-8895AA
           [8910., 8930.]]  # a skyline between 8910-8930AA
wav_cont = [[5560., 5565., 5590., 5595.],
            [6850., 6855., 6872.5, 6877.5],
            [8389., 8391., 8407.5, 8408.5],
            [8875., 8878., 8894., 8896.],
            [8910., 8912., 8930., 8935.]]

colname = ['wave0','e_wave0','lsigma','e_lsigma','vsigma','e_vsigma','flxcount','e_flxcount',
           'R','e_R','RMS','reduced_chi2']
df_skyline = pd.DataFrame(columns=colname)

n_fit = 1000

for i in tqdm.trange(np.shape(wav_fit)[0]):
    spx_fit = [np.abs(wav_new-wav_fit[i][0]).argmin(), np.abs(wav_new-wav_fit[i][1]).argmin()]

    spx_cont_left = [np.abs(wav_new-wav_cont[i][0]).argmin(), np.abs(wav_new-wav_cont[i][1]).argmin()]
    spx_cont_right = [np.abs(wav_new-wav_cont[i][2]).argmin(), np.abs(wav_new-wav_cont[i][3]).argmin()]

    cont_array = np.concatenate([d_sky[spx_cont_left[0]:spx_cont_left[1]],
                                 d_sky[spx_cont_right[0]:spx_cont_right[1]]], axis=0)
    cont = np.mean(cont_array)

    x_bin = wav_new[1]-wav_new[0]
    x_fit = wav_new[spx_fit[0]:spx_fit[1]+1]

    par0, par1, par2 = [], [], []
    for N in np.arange(n_fit):
        y_dat = np.random.normal(d_sky[spx_fit[0]:spx_fit[1]+1],
                                 np.sqrt(d_sky[spx_fit[0]:spx_fit[1]+1]))
        y_fit = y_dat-cont

        flx_scale0 = np.sum(y_dat-cont)*x_bin

        popt, pcov = curve_fit(gauss_cdf_scale, x_fit, y_fit,
                               bounds=([x_fit[0], x_bin, flx_scale0-100.0],
                                       [x_fit[-1], wav_fit[i][1]-wav_fit[i][0], flx_scale0+100.0]))
        par0.append(popt[0])
        par1.append(popt[1])
        par2.append(popt[2])

    lmu, e_lmu = np.mean(par0), np.std(par0)
    lsig, e_lsig = np.mean(par1), np.std(par1)
    vsig = c*lsig/lmu
    e_vsig = vsig * np.sqrt((e_lmu/lmu)**2.0 + (e_lsig/lsig)**2.0)
    flux, e_flux = np.mean(par2), np.std(par2)

    resol = lmu/(lsig*2.0*np.sqrt(2.0*np.log(2.0)))
    e_resol = resol * np.sqrt((e_lmu/lmu)**2.0 + (e_lsig/lsig)**2.0)

    y_var = d_sky[spx_fit[0]:spx_fit[1]+1]
    y_obs = d_sky[spx_fit[0]:spx_fit[1]+1]
    y_cal = cont + gauss_cdf_scale(x_fit, lmu, lsig, flux)

    rms = np.sqrt(np.sum((y_obs-y_cal)**2.0)/len(y_obs))
    rchisq = np.sum((y_obs-y_cal)**2.0 / y_var) / (len(y_obs)-3)

    df = pd.Series(data = {colname[0] : lmu,
                           colname[1] : e_lmu,
                           colname[2] : lsig,
                           colname[3] : e_lsig,
                           colname[4] : vsig,
                           colname[5] : e_vsig,
                           colname[6] : flux,
                           colname[7] : e_flux,
                           colname[8] : resol,
                           colname[9] : e_resol,
                           colname[10] : rms,
                           colname[11] : rchisq})

    df_skyline = df_skyline.append(df, ignore_index=True)


    ### Figure 2 : sky fitting
    fig1 = plt.figure(i+10, figsize=(8,4))
    ax1 = plt.subplot(1,1,1)
    ax1.set_position([0.10, 0.17, 0.85, 0.77])
    ax1.set_xlim([wav_fit[i][0]-25.0, wav_fit[i][1]+25.0])
    ax1.set_ylim([np.min(y_dat)-25.0, np.max(y_dat)+25.0])
    # ax1.set_xticks(np.linspace(5000, 10000, 6))
    # ax1.set_yticks(np.linspace(0, 1000, 6))
    # ax1.set_xticklabels(fontsize=13)
    # ax1.set_yticklabels(fontsize=13)
    ax1.set_xlabel(r'Wavelength $[\AA]$', fontsize=13)
    ax1.set_ylabel('Count', fontsize=13)
    ax1.tick_params(width=1.5, length=8.0)
    plt.minorticks_on()
    ax1.tick_params(width=1.5,length=5.0,which='minor')
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(1.5)

    ax1.plot(wav_new, d_sky, '-', color='tab:blue', linewidth=1.5, alpha=0.8)
    ax1.plot(wav_new, cont + gauss_cdf_scale(wav_new, lmu, lsig, flux),
             '-', color='tab:red', linewidth=2.0, alpha=0.6)
    ax1.axvline(x=lmu, linestyle='-', color='dimgray', linewidth=1.75, alpha=0.5)
    ax1.axvline(x=lmu-3.0*lsig, linestyle='--', color='dimgray', linewidth=1.75, alpha=0.5)
    ax1.axvline(x=lmu+3.0*lsig, linestyle='--', color='dimgray', linewidth=1.75, alpha=0.5)

    ax1.text(0.03, 0.95, r'$\lambda_{0}=$'+r'${0:.3f}\pm{1:.3f}\AA$'.format(lmu, e_lmu),
             ha='left', va='top', transform=ax1.transAxes, fontsize=10.0)
    ax1.text(0.03, 0.88, r'$\sigma_{v,inst}=$'+r'${0:.2f}\pm{1:.2f}~\rm km/s$'.format(vsig, e_vsig),
             ha='left', va='top', transform=ax1.transAxes, fontsize=10.0)
    ax1.text(0.03, 0.81, r'$R=$'+r'${0:.1f}\pm{1:.1f}$'.format(resol, e_resol),
             ha='left', va='top', transform=ax1.transAxes, fontsize=10.0)
    ax1.text(0.03, 0.74, 'RMS = '+'{0:.2e}'.format(rms),
             ha='left', va='top', transform=ax1.transAxes, fontsize=10.0)
    ax1.text(0.03, 0.67, r'$\chi_{v}^{2}=$'+'{0:.3f}'.format(rchisq),
             ha='left', va='top', transform=ax1.transAxes, fontsize=10.0)

    plt.savefig(dir_fig+'sky_fitting_{0:02d}.png'.format(i+1), dpi=300)
    plt.close()

df_skyline.to_pickle('df_skyline.pkl')


### Figure 3 : wavelength - resolution diagram
fig1 = plt.figure(50, figsize=(5,5))
ax1 = plt.subplot(1,1,1)
ax1.set_position([0.17, 0.17, 0.75, 0.75])
ax1.set_xticks(np.linspace(5000, 10000, 6))
ax1.set_yticks(np.linspace(400, 1600, 7))
ax1.set_xticklabels(['5000','6000','7000','8000','9000','10000'], fontsize=13)
ax1.set_yticklabels(['400','600','800','1000','1200','1400','1600'], fontsize=13)
ax1.set_xlim([4500.0, 9500.0])
ax1.set_ylim([500.0, 1400.0])
ax1.set_xlabel(r'Wavelength $[\AA]$ (observed frame)', fontsize=13)
ax1.set_ylabel('Spectral resolution', fontsize=13)
ax1.tick_params(width=1.5, length=8.0)
plt.minorticks_on()
ax1.tick_params(width=1.5,length=5.0,which='minor')
for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(1.5)

ax1.errorbar(df_skyline['wave0'].values, df_skyline['R'].values,
             xerr=df_skyline['e_wave0'].values, yerr=df_skyline['e_R'].values,
             fmt='o', ms=8.0, mew=1.5, mfc='deeppink', mec='black', ecolor='deeppink',
             capsize=0, capthick=2.0, elinewidth=2.0)

# Fitting 1
def temfunc(p, x):
    val = p[0]+p[1]*x
    return val

linear = Model(temfunc)

xfit = df_skyline['wave0'].values
yfit = df_skyline['R'].values
xefit = df_skyline['e_wave0'].values
yefit = df_skyline['e_R'].values

mydata = RealData(xfit, yfit, sx=xefit, sy=yefit)
myodr = ODR(mydata, linear, beta0=[0.1,1.0])
myoutput = myodr.run()
myoutput.pprint()

p = np.array([myoutput.beta[0], myoutput.beta[1]])
pe = np.array([myoutput.sd_beta[0], myoutput.sd_beta[1]])

ycal = temfunc(p, xfit)

ax1.plot(np.linspace(4500, 10000, num=11001, endpoint=True),
         temfunc(p, np.linspace(4500, 10000, num=11001, endpoint=True)),
         color='deeppink', linestyle='--', linewidth=1.75, alpha=0.8)

ax1.text(0.05, 0.92, r'$R=(%.1f\pm%.1f)+(%.2f\pm%.2f)\times\lambda_{\rm obs}$' %(p[0], pe[0], p[1], pe[1]),
         ha='left', va='top', transform=ax1.transAxes, fontsize=10.0, color='crimson')
ax1.text(0.05, 0.85, 'RMS = {0:.2f}'.format(np.sqrt(np.sum((yfit-ycal)**2.0)/len(yfit))),
         ha='left', va='top', transform=ax1.transAxes, fontsize=10.0, color='crimson')
# func_res = interpolate.interp1d(x_inter, y_inter, kind='linear', fill_value='extrapolate')

# ax1.plot(np.linspace(4500, 10000, num=11001, endpoint=True), func_res(np.linspace(4500, 10000, 5501)),
#          color='deeppink', linestyle='--', linewidth=1.75, alpha=0.7)

plt.savefig(dir_fig+'sky_resolution.png', dpi=300)
plt.close()


# ----- Printing resolution values ----- #
elines = ['[OII]3727+3729', 'Hbeta' , '[OIII]4959', '[OIII]5007',
          '[NII]6548', 'Halpha', '[NII]6584', '[SII]6717+6731']
ewav = [(3727.092+3729.875)/2, 4862.680, 4960.295, 5008.240,
        6549.860, 6564.610, 6585.270, (6718.29+6732.67)/2]

print('\n')

colname = ['elines', 'ewav', 'R', 'e_R', 'vd_inst', 'e_vd_inst']
df_resol = pd.DataFrame(columns=colname)

for i in np.arange(len(elines)):
    print('# ----- At '+elines[i]+' line ----- #')
    eR = temfunc(p, ewav[i]*(1.0+ic.redshift))
    e_eR = np.sqrt(pe[0]**2.0 + (ewav[i]*(1.0+ic.redshift)*pe[1])**2.0)
    print('R : {0:.1f} +/- {1:.1f}'.format(eR, e_eR))
    vd_inst = c / (2.0*np.sqrt(2.0*np.log(2.0))*eR)
    e_vd_inst = c*e_eR / (2.0*np.sqrt(2.0*np.log(2.0))*eR*eR)
    # ld_inst = ewav[i]*(1.0+ic.redshift) / (2.0*np.sqrt(2.0*np.log(2.0))*eR)
    # vd_inst = c*ld_inst/(ewav[i]*(1.0+ic.redshift))
    print('vd_inst : {0:.2f} +/- {1:.2f} km/s'.format(vd_inst, e_vd_inst))
    print('\n')

    df = pd.Series(data = {colname[0] : elines[i],
                           colname[1] : ewav[i],
                           colname[2] : eR,
                           colname[3] : e_eR,
                           colname[4] : vd_inst,
                           colname[5] : e_vd_inst})
    df_resol = df_resol.append(df, ignore_index=True)

df_resol.to_pickle('df_resol.pkl')


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))

