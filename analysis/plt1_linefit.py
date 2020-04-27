#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 9 14:45:23 2020

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import glob, os
import init_cfg as ic
from astropy.io import fits
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import tqdm
from scipy.special import erf
import pandas as pd


def gauss_cdf_scale(x, mu, sigma, flux_scale):
    dx = x[1] - x[0]
    v1 = erf((x-mu+0.5*dx)/(np.sqrt(2.0)*sigma))
    v2 = erf((x-mu-0.5*dx)/(np.sqrt(2.0)*sigma))
    return flux_scale*(v1-v2)/(2.0*dx)

def multi2_gauss_cdf_scale(x, *pars):
    g1 = gauss_cdf_scale(x, pars[0], pars[1], pars[2])
    g2 = gauss_cdf_scale(x, pars[3], pars[4], pars[5])
    return g1+g2

def multi3_gauss_cdf_scale(x, *pars):
    g1 = gauss_cdf_scale(x, pars[0], pars[1], pars[2])
    g2 = gauss_cdf_scale(x, pars[3], pars[4], pars[5])
    g3 = gauss_cdf_scale(x, pars[6], pars[7], pars[8])
    return g1+g2+g3


# ----- Directory ----- #
dir_fig = '/data/jlee/DATA/Gemini/Programs/GN-2019A-Q-215/analysis/diagram/linefits/'

dir_lines = '/data/jlee/DATA/Gemini/Programs/GN-2019A-Q-215/analysis/lines/'
glob_lines = glob.glob(dir_lines+'*')

emi_lines = []
for i in np.arange(len(glob_lines)):
	emi_lines.append(glob_lines[i].split(dir_lines)[1])
emi_lines = sorted(emi_lines)


for i in np.arange(len(emi_lines)):
    exec('dir_'+emi_lines[i]+' = "'+dir_lines+emi_lines[i]+'/"')

# emi_lines
# ['Halpha',
#  'Hbeta',
#  'NII6548',
#  'NII6584',
#  'OII3727',
#  'OIII4959',
#  'OIII5007',
#  'SII6717',
#  'SII6731']


# ----- Basic parameters ----- #
redshift = 0.353
dist_lum = 1875.5e+6    # pc
c = 2.99792e+5    # km/s
# Angstrom (SDSS)
wav_lines = [6564.61, 4862.68, 6549.86, 6585.27,
             0.5 * (3727.092 + 3729.875),
             4960.295, 5008.240, 6718.29, 6732.67]

for i in np.arange(len(emi_lines)):
    exec('wav_'+emi_lines[i]+' = '+str(wav_lines[i]))


# ----- Reading 2D data ----- #
for l in np.arange(len(emi_lines)):
    exec("dat_dir = dir_"+emi_lines[l])
    dat2D_list = sorted(glob.glob(dat_dir+'*_2D.fits'))
    for i in np.arange(len(dat2D_list)):
        dat_2D = emi_lines[l]+'_'+dat2D_list[i].split('/')[-1].split('.fits')[0]
        exec(dat_2D+" = fits.getdata('"+dat2D_list[i]+"', ext=1)")
        if (dat_2D in [emi_lines[l]+'_'+'lmu_2D', emi_lines[l]+'_'+'lsig_2D',
                       emi_lines[l]+'_'+'vsig_2D', emi_lines[l]+'_'+'flx_2D']):
            exec('e_'+dat_2D+" = fits.getdata('"+dat2D_list[i]+"', ext=2)")


# ----- Reading the cube ----- #
fin_cb = 'sfcube_3D.fits'

hd0 = fits.getheader(fin_cb, ext=0)
d_sci2, h_sci = fits.getdata(fin_cb, ext=1, header=True)
d_var2, h_var = fits.getdata(fin_cb, ext=2, header=True)

wav = np.linspace(start=h_sci['CRVAL3'],
                  stop=h_sci['CRVAL3']+(h_sci['NAXIS3']-1)*h_sci['CD3_3'],
                  num=h_sci['NAXIS3'], endpoint=True)

wav_rest = wav/(1.0+redshift)
d_sci2 *= (1.0+redshift)
d_var2 *= (1.0+redshift)**2.0


# ----- Plotting figures ----- #
fit_2D = ['cont_2D', 'lmu_2D', 'rchisq_2D', 'snr_2D',
          'flx_2D', 'lsig_2D', 'rms_2D', 'vsig_2D']

for l in np.arange(len(emi_lines)):

	for j in np.arange(len(fit_2D)):
		var = emi_lines[l]+'_'+fit_2D[j]

		fig, ax = plt.subplots()
		ax.set_xlim([0.5,68.5])
		ax.set_ylim([0.5,49.5])
		ax.set_xlabel('pixel')
		ax.set_ylabel('pixel')

		exec('v_low, v_high = np.percentile('+var+', [2.0, 98.0])')
		exec("im = ax.imshow("+var+", cmap='gray_r', vmin=v_low, vmax=v_high, origin='lower', extent=[0.5,68.5,0.5,49.5])")
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		plt.colorbar(im, cax=cax)

		plt.savefig(dir_fig+var+'.pdf')
		plt.close()


# ----- Selection criteria ----- #
df_resol = pd.read_pickle('df_resol.pkl')
vd_inst = df_resol.loc[df_resol['elines'] == 'Halpha']['vd_inst'].values[0]

crit_2D = ((Halpha_flx_2D > 0.0) & \
           (Halpha_lmu_2D > 6560.0) & \
           (Halpha_lmu_2D < 6570.0) & \
           (Halpha_vsig_2D > vd_inst) & \
           (Halpha_snr_2D > 3.0))
crit_2D[:,:20-1] = False


# ----- H alpha flux map ----- #
plt_Data = Halpha_flx_2D*crit_2D

fig, ax = plt.subplots()
ax.set_xlim([0.5,68.5])
ax.set_ylim([0.5,49.5])
ax.set_xlabel('pixel') 
ax.set_ylabel('pixel')

v_low, v_high = np.percentile(plt_Data, [2.0, 98.0])
im = ax.imshow(plt_Data, cmap='gray_r', vmin=v_low, vmax=v_high, origin='lower', extent=[0.5,68.5,0.5,49.5])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.savefig(dir_fig+'Halpha_fluxmap.pdf')
plt.close()


# ----- SFR (H alpha) map ----- #
EBV_gal = 0.026
EBV_int_2D = 1.97*np.log10((Halpha_flx_2D/Hbeta_flx_2D)/2.86)
k_Halpha = 3.32
A_Halpha = k_Halpha * (EBV_int_2D + EBV_gal)

L_Halpha = 1.0e-15 * Halpha_flx_2D * 10.0**(0.4*A_Halpha) * (4.0*np.pi*(dist_lum*3.086e+18)**2.0)
plt_Data = L_Halpha * 4.6e-42 * crit_2D

fig, ax = plt.subplots()
ax.set_xlim([0.5,68.5])
ax.set_ylim([0.5,49.5])
ax.set_xlabel('pixel') 
ax.set_ylabel('pixel')

v_low, v_high = np.percentile(plt_Data, [2.0, 98.0])
im = ax.imshow(plt_Data, cmap='gray_r', vmin=v_low, vmax=v_high, origin='lower', extent=[0.5,68.5,0.5,49.5])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.savefig(dir_fig+'Halpha_SFR.pdf')
plt.close()


# ----- H alpha radial velocity distribution map ----- #
ymax_idx, xmax_idx = np.unravel_index(Halpha_flx_2D.argmax(), Halpha_flx_2D.shape)

# plt_Data = (c*(Halpha_lmu_2D-wav_Halpha)/wav_Halpha)*crit_2D
plt_Data = (c*(Halpha_lmu_2D-Halpha_lmu_2D[ymax_idx, xmax_idx])/Halpha_lmu_2D[ymax_idx, xmax_idx])*crit_2D
v_low, v_high = np.percentile(plt_Data[crit_2D*1 == 1], [2.0, 98.0])

plt_Data[crit_2D*1 == 0] = np.nan

fig, ax = plt.subplots()
ax.set_xlim([0.5,68.5])
ax.set_ylim([0.5,49.5])
ax.set_xlabel('pixel') 
ax.set_ylabel('pixel')

im = ax.imshow(plt_Data, cmap='rainbow', vmin=v_low, vmax=v_high, origin='lower', extent=[0.5,68.5,0.5,49.5])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.savefig(dir_fig+'Halpha_rvmap.pdf')
plt.close()


# ----- H alpha velocity dispersion map ----- #
plt_Data = np.sqrt(Halpha_vsig_2D**2.0 - vd_inst**2.0)*crit_2D
# plt_Data[np.isnan(plt_Data) == True] = 0.0
v_low, v_high = np.percentile(plt_Data[crit_2D*1 == 1], [2.0, 98.0])

plt_Data[crit_2D*1 == 0] = np.nan

fig, ax = plt.subplots()
ax.set_xlim([0.5,68.5])
ax.set_ylim([0.5,49.5])
ax.set_xlabel('pixel') 
ax.set_ylabel('pixel')

im = ax.imshow(plt_Data, cmap='rainbow', vmin=v_low, vmax=v_high, origin='lower', extent=[0.5,68.5,0.5,49.5])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.savefig(dir_fig+'Halpha_vdisp.pdf')
plt.close()


# ----- H alpha / H beta ratio map ----- #
plt_Data = (Halpha_flx_2D/Hbeta_flx_2D)*crit_2D
v_low, v_high = np.percentile(plt_Data[crit_2D*1 == 1], [2.0, 98.0])

plt_Data[crit_2D*1 == 0] = np.nan

fig, ax = plt.subplots()
ax.set_xlim([0.5,68.5])
ax.set_ylim([0.5,49.5])
ax.set_xlabel('pixel') 
ax.set_ylabel('pixel')

im = ax.imshow(plt_Data, cmap='rainbow', vmin=v_low, vmax=v_high, origin='lower', extent=[0.5,68.5,0.5,49.5])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.savefig(dir_fig+'Hab_ratio.pdf')
plt.close()


# ----- N2 / H alpha ratio map ----- #
plt_Data = (NII6584_flx_2D/Halpha_flx_2D)*crit_2D
v_low, v_high = np.percentile(plt_Data[crit_2D*1 == 1], [2.0, 98.0])

plt_Data[crit_2D*1 == 0] = np.nan

fig, ax = plt.subplots()
ax.set_xlim([0.5,68.5])
ax.set_ylim([0.5,49.5])
ax.set_xlabel('pixel') 
ax.set_ylabel('pixel')

im = ax.imshow(plt_Data, cmap='rainbow', vmin=v_low, vmax=v_high, origin='lower', extent=[0.5,68.5,0.5,49.5])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.savefig(dir_fig+'N2Ha_ratio.pdf')
plt.close()


# ----- O3 / H beta ratio map ----- #
plt_Data = (OIII5007_flx_2D/Hbeta_flx_2D)*crit_2D
v_low, v_high = np.percentile(plt_Data[crit_2D*1 == 1], [2.0, 98.0])

plt_Data[crit_2D*1 == 0] = np.nan

fig, ax = plt.subplots()
ax.set_xlim([0.5,68.5])
ax.set_ylim([0.5,49.5])
ax.set_xlabel('pixel') 
ax.set_ylabel('pixel')

im = ax.imshow(plt_Data, cmap='rainbow', vmin=v_low, vmax=v_high, origin='lower', extent=[0.5,68.5,0.5,49.5])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.savefig(dir_fig+'O3Hb_ratio.pdf')
plt.close()


# ----- S2 / H alpha ratio map ----- #
plt_Data = (SII6717_flx_2D/Halpha_flx_2D)*crit_2D
v_low, v_high = np.percentile(plt_Data[crit_2D*1 == 1], [2.0, 98.0])

plt_Data[crit_2D*1 == 0] = np.nan

fig, ax = plt.subplots()
ax.set_xlim([0.5,68.5])
ax.set_ylim([0.5,49.5])
ax.set_xlabel('pixel') 
ax.set_ylabel('pixel')

im = ax.imshow(plt_Data, cmap='rainbow', vmin=v_low, vmax=v_high, origin='lower', extent=[0.5,68.5,0.5,49.5])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.savefig(dir_fig+'S2Ha_ratio.pdf')
plt.close()


# BPT diagram
fig, ax = plt.subplots()
ax.set_xlim([-1.4, 0.5])
ax.set_ylim([-1.4, 1.5])
ax.set_xlabel(r'log([NII]$\lambda6584$/H$\alpha$)') 
ax.set_ylabel(r'log([OIII]$\lambda5007$/H$\beta$)')

x_Dat = np.log10(((NII6584_flx_2D/Halpha_flx_2D)*crit_2D))
y_Dat = np.log10(((OIII5007_flx_2D/Hbeta_flx_2D)*crit_2D))

BPT_SFG = ((crit_2D*1 == 1) & (y_Dat < 1.3+0.61/x_Dat))
BPT_comp = ((crit_2D*1 == 1) & (y_Dat > 1.3+0.61/x_Dat) & (y_Dat < 1.19+0.61/(x_Dat-0.47)))
BPT_AGN = ((crit_2D*1 == 1) & (y_Dat > 1.19+0.61/(x_Dat-0.47)) & \
           (y_Dat-1.4899 > ((0.2949- -0.2135)/(1.4899-0.4548))*(x_Dat-0.2949)))
BPT_LINER = ((crit_2D*1 == 1) & (y_Dat > 1.19+0.61/(x_Dat-0.47)) & \
             (y_Dat-1.4899 < ((0.2949- -0.2135)/(1.4899-0.4548))*(x_Dat-0.2949)))

x_Dat_flat = x_Dat.flatten()
y_Dat_flat = y_Dat.flatten()

ax.plot(x_Dat_flat[BPT_SFG.flatten()], y_Dat_flat[BPT_SFG.flatten()], 'o', ms=2.75, color='dodgerblue', alpha=0.8)
ax.plot(x_Dat_flat[BPT_comp.flatten()], y_Dat_flat[BPT_comp.flatten()], 'o', ms=2.75, color='darkorange', alpha=0.8)
ax.plot(x_Dat_flat[BPT_AGN.flatten()], y_Dat_flat[BPT_AGN.flatten()], 'o', ms=2.75, color='red', alpha=0.8)
ax.plot(x_Dat_flat[BPT_LINER.flatten()], y_Dat_flat[BPT_LINER.flatten()], 'o', ms=2.75, color='green', alpha=0.8)


# Kauffmann+03
xx = np.linspace(-1.5, -0.01, 1000)
ax.plot(xx, 1.3+0.61/xx, 'k-', linewidth=1.25, alpha=0.75)

# Kewley+01
xx = np.linspace(-1.5, 0.45, 1000)
ax.plot(xx, 1.19+0.61/(xx-0.47), 'k:', linewidth=1.25, alpha=0.75)

# AGN-LINER boundary (digitized from GASP XV. paper)
ax.plot([-0.2135, 0.4548], [0.2949, 1.4899], 'k--', linewidth=1.25, alpha=0.75)


plt.savefig(dir_fig+'BPT.pdf')
plt.close()


# ----- BPT spatial map ----- #
# BPT_fld = np.empty((d_sci2.shape[1], d_sci2.shape[2]))
# BPT_fld[:] = np.nan

# BPT_fld[BPT_SFG] = 
# BPT_fld[BPT_comp] = 
# BPT_fld[BPT_AGN] = 
# BPT_fld[BPT_LINER] = 



# plt_Data = (SII6717_flx_2D/Halpha_flx_2D)*crit_2D
# v_low, v_high = np.percentile(plt_Data[crit_2D*1 == 1], [2.0, 98.0])

# plt_Data[crit_2D*1 == 0] = np.nan

fig, ax = plt.subplots()
ax.set_xlim([0.5,68.5])
ax.set_ylim([0.5,49.5])
ax.set_xlabel('pixel') 
ax.set_ylabel('pixel')

for i in np.arange(np.sum(BPT_SFG == True)):
	rect = plt.Rectangle((np.argwhere(BPT_SFG == True)[i,1], np.argwhere(BPT_SFG == True)[i,0]),
		                 1.0, 1.0, facecolor='dodgerblue', alpha=0.8)
	ax.add_patch(rect)

for i in np.arange(np.sum(BPT_comp == True)):
	rect = plt.Rectangle((np.argwhere(BPT_comp == True)[i,1], np.argwhere(BPT_comp == True)[i,0]),
		                 1.0, 1.0, facecolor='darkorange', alpha=0.8)
	ax.add_patch(rect)

for i in np.arange(np.sum(BPT_AGN == True)):
	rect = plt.Rectangle((np.argwhere(BPT_AGN == True)[i,1], np.argwhere(BPT_AGN == True)[i,0]),
		                 1.0, 1.0, facecolor='red', alpha=0.8)
	ax.add_patch(rect)

for i in np.arange(np.sum(BPT_LINER == True)):
	rect = plt.Rectangle((np.argwhere(BPT_LINER == True)[i,1], np.argwhere(BPT_LINER == True)[i,0]),
		                 1.0, 1.0, facecolor='green', alpha=0.8)
	ax.add_patch(rect)


plt.savefig(dir_fig+'BPTmap.pdf')
plt.close()



# Printing the running time
print('\n')
print('--- %s seconds ---' %(time.time()-start_time))