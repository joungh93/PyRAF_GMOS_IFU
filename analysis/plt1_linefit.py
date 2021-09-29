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
from astropy.table import Table 
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
dir_fig = ic.cpath+"diagram/linefits/"
diG = ic.dir_redux
dir_lines = ic.cpath+"lines3/"
glob_lines = glob.glob(dir_lines+'*')

emi_lines = []
for i in np.arange(len(glob_lines)):
    emi_lines.append(glob_lines[i].split(dir_lines)[1])
emi_lines = sorted(emi_lines)
emi_lines.remove('check')


for i in np.arange(len(emi_lines)):
    exec('dir_'+emi_lines[i]+' = "'+dir_lines+emi_lines[i]+'/"')

# emi_lines
# ['Halpha', 'Hbeta', 'NII6548', 'NII6584',
#  'OII3727', 'OIII4959', 'OIII5007',
#  'SII6717', 'SII6731']

name_elines = [r"${\rm H\alpha}$", r"${\rm H\beta}$", r"${\rm [NII]\lambda6548}$", r"${\rm [NII]\lambda6584}$",
               r"${\rm [OII]\lambda\lambda3727,3729}$", r"${\rm [OIII]\lambda4959}$", r"${\rm [OIII]\lambda5007}$",
               r"${\rm [SII]\lambda6717}$", r"${\rm [SII]\lambda6731}$"]


# ----- Basic parameters ----- #
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
redshift = ic.redshift
dist_lum = cosmo.luminosity_distance(redshift).value*1.0e+6    # pc
c = 2.99792e+5    # km/s
ang_scale = cosmo.kpc_proper_per_arcmin(redshift).value / 60.    # kpc/arcsec
pixel_scale = ic.pixel_scale    # arcsec/pix

# Angstrom (SDSS)
wav_lines = [6564.61, 4862.68, 6549.86, 6585.27,
             3727.092, 4960.295, 5008.240,
             6718.29, 6732.67]
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
hdr1 = fits.getheader(ic.cube_list[0], ext=0)
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
X_coord = pixel_scale*(np.arange(Halpha_flux_2D.shape[1], step=1)-Halpha_flux_2D.shape[1]/2)
Y_coord = pixel_scale*(np.arange(Halpha_flux_2D.shape[0], step=1)-Halpha_flux_2D.shape[0]/2)

flx_Data = copy.deepcopy(Halpha_flux_2D)
snr_cnd = ((Halpha_snr_2D < 3.0) | (Halpha_snrpix_2D < 5.0))
# sig_cnd = (Halpha_sigma_2D < lsig_llim)
rchi25, rchi50, rchi75 = np.percentile(Halpha_rchisq_2D[Halpha_rchisq_2D > 0.],
                                       [25.0, 50.0, 75.0])
rchisq_cnd = (Halpha_rchisq_2D > 50.)
flx25, flx50, flx75 = np.percentile(Halpha_flux_2D[Halpha_flux_2D > 0.],
                                    [25.0, 50.0, 75.0])
flx_cnd = (Halpha_flux_2D < flx25)
zero_cnd = (snr_cnd | rchisq_cnd)

flx_Data[zero_cnd] = 0.
sflx = ndimage.gaussian_filter(flx_Data, sigma=(1.25,1.25), order=0)

sig = np.std(sflx)
lvs = [0.25*sig, 0.5*sig, 1.*sig, 3.*sig, 5.*sig, 7.5*sig]
lws = tuple(np.repeat(2.5, len(lvs)-1))
cs = tuple(['gray']*(len(lvs)-1))


# ----- Custom functions ----- #
pltFlags = {'x0':3.00, 'y0':1.75, 'sign':1, 'L':0.6, 'theta0':gpa*(np.pi/180.0),
            'xN':2.00, 'yN':2.10, 'xE':2.40, 'yE':0.80, 'legend_position':'lower left'}


def plot_2Dmap(plt_Data, title, v_low, v_high, out, cmap='gray_r',
               add_cb=True, add_ct=True, add_or=True, add_sc=True,
               cb_label=None, add_legend=True, legend_position='lower left',
               x0=-2.75, y0=1.25, sign=-1, L=0.6, theta0=gpa*(np.pi/180.0),
               xN=-1.90, yN=1.25, xE=-2.95, yE=2.10):

    fig, ax = plt.subplots(1, 1, figsize=(8,5))
    plt.suptitle(title, x=0.5, ha='center', y=0.96, va='top',
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

    im = ax.imshow(plt_Data, cmap=cmap, vmin=v_low, vmax=v_high,
                   aspect='equal', extent=[-3.4,3.4,-2.45,2.45])

    if add_cb:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = plt.colorbar(im, cax=cax)
        cb.set_label(cb_label, size=12.0, labelpad=15.0)
        cb.ax.tick_params(labelsize=12.0)

    if add_ct:
        ax.contour(X_coord, Y_coord[::-1], sflx, levels=lvs, linewidths=lws, colors=cs, alpha=0.6)
        if add_legend:
            p0, = ax.plot(-100.0, -100.0, '-', linewidth=2.5, color='gray', alpha=0.6,
                          label=r"H${\rm \alpha}$ flux contour")
            ax.legend(handles=[p0], fontsize=13.0, loc=legend_position,
                      handlelength=2.5, frameon=True, borderpad=0.8,
                      framealpha=0.8, edgecolor='gray')

    if add_or:
        ax.arrow(x0+sign*0.025, y0, L*np.sin(theta0), L*np.cos(theta0), width=0.06,
                 head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
        ax.arrow(x0, y0+sign*0.025, -L*np.cos(theta0), L*np.sin(theta0), width=0.06,
                 head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
        ax.text(xN, yN, 'N', fontsize=15.0, fontweight='bold', color='blueviolet')
        ax.text(xE, yE, 'E', fontsize=15.0, fontweight='bold', color='blueviolet')
    
    if add_sc:
        kpc5 = 5.0 / ang_scale
        ax.arrow(2.0, -1.85, kpc5, 0., width=0.07, head_width=0., head_length=0.,
                  fc='blueviolet', ec='blueviolet', alpha=0.9)
        ax.text(2.1, -2.2, '5 kpc', fontsize=15.0, fontweight='bold', color='blueviolet')

    plt.savefig(out+'.pdf', dpi=300)
    plt.savefig(out+'.png', dpi=300)
    plt.close()


def weighted_mean(data, weights, e_data=None, e_weights=None):
    wm = np.average(data, weights=weights)

    if ((e_data is not None) & (e_weights is not None)):
        wht_sum = np.sum(weights)
        e_wht_sum = np.sqrt(np.sum(e_weights**2.))

        Aj = weights * data
        e_Aj = Aj * np.sqrt((e_weights/weights)**2 + (e_data/data)**2)

        Cj = Aj / wht_sum
        e_Cj = Cj * np.sqrt((e_Aj/Aj)**2 + (e_wht_sum/wht_sum)**2)

        e_wm = np.sum(e_Cj)

        res = [wm, e_wm]

    else:
        res = wm
    
    return res


def get_line_ratio(flux1, flux2, snr1, snr2, e_flux1=None, e_flux2=None):
    # flux1, flux2 : flux data 2D array
    line_ratio = flux1 / flux2
    snr2_cnd = ((snr1 < 3.0) | (snr2 < 3.0))
    zero2_cnd = (zero_cnd | snr2_cnd)
    line_ratio[zero2_cnd] = 0.
    line_ratio[line_ratio == 0.] = np.nan
    line_ratio[np.isinf(line_ratio) == True] = np.nan
    if ((e_flux1 is not None) & (e_flux2 is not None)):
        e_line_ratio = line_ratio * np.sqrt((e_flux1/flux1)**2 + (e_flux2/flux2)**2)
    else:
        e_line_ratio = np.zeros_like(line_ratio)
    return [line_ratio, e_line_ratio]
# ---------------------------- #


# ----- START: Flux maps ----- #
for l in np.arange(len(emi_lines)):
    exec("snr2_cnd = ("+emi_lines[l]+"_snr_2D < 3.0)")
    zero2_cnd = (zero_cnd | snr2_cnd)

    exec(emi_lines[l]+"_flux_2D[zero2_cnd] = 0.")
    exec("e_"+emi_lines[l]+"_flux_2D[zero2_cnd] = 0.")
    exec("plt_Data = "+emi_lines[l]+"_flux_2D")
    v_low, v_high = np.percentile(plt_Data[plt_Data > 0.], [10.0, 90.0])
    plot_2Dmap(plt_Data, name_elines[l]+" flux map", v_low, v_high,
               dir_fig+"Map_flux_"+emi_lines[l],
               cb_label=r'Flux [${\rm 10^{-15}~erg~s^{-1}~cm^{-2}~\AA^{-1}}$]', **pltFlags)
# ----- END: Flux maps ----- #


# ----- START: S/N maps ----- #
for l in np.arange(len(emi_lines)):
    exec("plt_Data = "+emi_lines[l]+"_snrpix_2D")
    plt_Data[zero_cnd] = 0.
    v_low, v_high = np.percentile(plt_Data[plt_Data > 0.], [10.0, 90.0])

    plot_2Dmap(plt_Data, name_elines[l]+" S/N map (per pixel)", 0., v_high,
               dir_fig+"Map_snr_"+emi_lines[l],
               cb_label="Signal-to-noise ratio per pixel", **pltFlags)
# ----- END: S/N maps ----- #


# ----- START: Radial velocity distribution (H alpha) map ----- #
ymax_idx, xmax_idx = np.unravel_index(Halpha_flux_2D.argmax(), Halpha_flux_2D.shape)
Halpha_mu0 = Halpha_mu_2D[ymax_idx, xmax_idx]

plt_Data = (c*(Halpha_mu_2D-Halpha_mu0)/Halpha_mu0)
plt_Data[zero_cnd] = np.nan
Rvd = plt_Data

v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
plot_2Dmap(plt_Data, "Radial velocity map", v_low-25.0, v_high+25.0,
           dir_fig+"Map_rv_Halpha", cmap='rainbow',
           cb_label=r"Relative velocity [${\rm km~s^{-1}}$]", **pltFlags)
# ----- END: Radial velocity distribution (H alpha) map ----- #


# ----- START: Velocity dispersion map ----- #
emi_lines2 = ["Halpha", "Hbeta", "OII3727", "OIII5007", "SII6717"]
name_elines2 = [r"${\rm H\alpha}$", r"${\rm H\beta}$",
                r"${\rm [OII]\lambda\lambda3727,3729}$", r"${\rm [OIII]\lambda5007}$",
                r"${\rm [SII]\lambda6717}$"]
wav_lines2 = [6564.61, 4862.68, 3727.092, 5008.240, 6718.29]
for i, l in enumerate(emi_lines2):
    Rmax1 = par[0] + par[1]*wav_lines2[i]*(1+redshift)
    vsig1 = c / (2.0*np.sqrt(2.0*np.log(2.0))*Rmax1)
    exec("vsig2 = "+l+"_vsig_2D")
    exec("e_vsig2 = "+l+"_e_vsig_2D")
    plt_Data = np.sqrt(vsig2**2.0 - vsig1**2.0)
    e_plt_Data = np.abs(vsig2 / plt_Data) * e_vsig2
    plt_Data[((vsig2 > 0.) & (vsig2 <= vsig1))] = 0.
    plt_Data[zero_cnd] = np.nan
    e_plt_Data[zero_cnd] = np.nan
    Vdd, e_Vdd = plt_Data, e_plt_Data

    v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])

    plot_2Dmap(plt_Data, r"Velocity dispersion map ("+name_elines2[i]+")",
               v_low-25.0, v_high+25.0, dir_fig+"Map_vd_"+l, cmap='rainbow',
               cb_label=r"Velocity dispersion [${\rm km~s^{-1}}$]", **pltFlags)

    # Flux-weighted mean
    val_vdd = ((np.isnan(plt_Data) == False) & (plt_Data > 0.))
    fwm_vdisp, e_fwm_vdisp = weighted_mean(data=plt_Data[val_vdd], weights=Halpha_flux_2D[val_vdd],
                                           e_data=e_plt_Data[val_vdd], e_weights=e_Halpha_flux_2D[val_vdd])

    print("Flux-weighted mean of velocity dispersion ("+l+f") : {fwm_vdisp:.2f} +/- {e_fwm_vdisp:.2f} km/s")

    # Error map
    plot_2Dmap(plt_Data/e_plt_Data, r"Velocity dispersion S/N map ("+name_elines2[i]+")", 0.0, 3.0,
               dir_fig+"Map_evd_"+l, cmap='rainbow', **pltFlags)
# ----- END: Velocity dispersion map ----- #


# ----- START: H alpha / H beta flux ratio map ----- #
plt_Data, e_plt_Data = get_line_ratio(Halpha_flux_2D, Hbeta_flux_2D, Halpha_snr_2D, Hbeta_snr_2D,
                                      e_Halpha_flux_2D, e_Hbeta_flux_2D)
v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
plot_2Dmap(plt_Data, r"${\rm H\alpha/H\beta}$ flux ratio map", 2.86, 1.10*v_high,
           dir_fig+"Line_ratio_Hab", cmap='rainbow', **pltFlags)

# Flux-weighted mean
val_Hab = (np.isnan(plt_Data) == False)
fwm_Hab, e_fwm_Hab = weighted_mean(data=plt_Data[val_Hab], weights=Halpha_flux_2D[val_Hab],
                                   e_data=e_plt_Data[val_Hab], e_weights=e_Halpha_flux_2D[val_Hab])
print(f"Flux-weighted mean of Ha/Hb flux ratio : {fwm_Hab:.3f} +/- {e_fwm_Hab:.3f}")

# Error map
plot_2Dmap(plt_Data/e_plt_Data, r"${\rm H\alpha/H\beta}$ flux ratio S/N map", 0.0, 3.0,
           dir_fig+"Line_ratio_eHab", cmap='rainbow', **pltFlags)
# ----- END: H alpha / H beta flux ratio map ----- #


# ----- START: SFR (H alpha) map ----- #
def compute_SFR(Ha_flux, Hab_ratio, luminosity_distance,
                e_Ha_flux=0.0, e_Hab_ratio=0.0, EBV_gal=0.0, apply_C00=True):
    '''
    luminosity_distance: [in pc]
    '''

    # Please refer to test_dustlaws.ipynb.
    if apply_C00:
        k_Ha, k_Hb, k_V = 3.3248, 4.5965, 4.0522
    else:  # Cardelli+89 extinction law
        k_Ha, k_Hb, k_V = 2.5342, 3.6076, 3.1000

    EBV_int = (-2.5 / (k_Ha-k_Hb)) * np.log10(Hab_ratio / 2.86)
    EBV = EBV_int + EBV_gal
    A_Ha = k_Ha * EBV
    e_A_Ha = k_Ha * (-2.5 / (k_Ha-k_Hb)) * e_Hab_ratio/(Hab_ratio*np.log(10.0))
    A_V = (k_V / k_Ha) * A_Ha

    L_Ha = 1.0e-15 * Ha_flux * 10.0**(0.4*A_Ha) * (4.0*np.pi*(luminosity_distance*3.086e+18)**2.0)
    e_L_Ha = L_Ha * np.sqrt((e_Ha_flux/Ha_flux)**2. + (e_A_Ha/A_Ha)**2.)

    SFR = L_Ha * 4.6e-42
    e_SFR = e_L_Ha * 4.6e-42
    e_SFR[SFR == 0.] = 0.

    return [SFR, e_SFR, A_V]

Hab = np.ones_like(plt_Data) * fwm_Hab    # Assumption (to be revised later)
EBV_gal = 0.012
SFR, e_SFR, A_V = compute_SFR(Halpha_flux_2D, Hab, dist_lum, 
                              e_Halpha_flux_2D, e_fwm_Hab, EBV_gal=EBV_gal, apply_C00=False)
val_SFR = (SFR > 0.)
print(f"SFR sum : {np.sum(SFR[val_SFR]):.2f} +/- {np.sqrt(np.sum(e_SFR[val_SFR]**2)):.2f} Mo/yr")
print(f"Mean V-magnitude extinction : {np.mean(A_V[val_SFR]):.3f} mag")

plt_Data = SFR / (pixel_scale*ang_scale)**2.
SFRD = plt_Data

v_low, v_high = np.percentile(plt_Data, [1.0, 99.0])
plot_2Dmap(plt_Data, "SFR map", v_low, v_high,
           dir_fig+"Map_SFR_Halpha", cmap='gray_r',
           cb_label=r"SFR density [$M_{\odot}~{\rm yr^{-1}~kpc^{-2}}$]", **pltFlags)
# ----- END: SFR (H alpha) map ----- #


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
    phot_table = apphot(data=SFR, apertures=ap, error=e_SFR, mask=np.logical_not(val_SFR))

    return [phot_table['aperture_sum'].data[0], phot_table['aperture_sum_err'].data[0]]

# # 1. HST
# SFR_disk_hst, e_SFR_disk_hst = cal_sfr_disk("HST_boundary_1sig_transformed.reg")
# print(f"SFR disk (HST) : {SFR_disk_hst:.2f} +/- {e_SFR_disk_hst:.2f} Mo/yr")

# 2. Gemini
SFR_disk_gem, e_SFR_disk_gem = cal_sfr_disk("GMOS_boundary_1sig.reg")
print(f"SFR disk (Gemini) : {SFR_disk_gem:.2f} +/- {e_SFR_disk_gem:.2f} Mo/yr")
# ----- END: calculating tail SFR ----- #


# ----- START: [NII]6584 / H alpha flux ratio map ----- #
plt_Data, e_plt_Data = get_line_ratio(NII6584_flux_2D, Halpha_flux_2D, NII6584_snr_2D, Halpha_snr_2D,
                                      e_NII6584_flux_2D, e_Halpha_flux_2D)
v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
plot_2Dmap(plt_Data, r"${\rm [NII]\lambda 6584/H\alpha}$ flux ratio map", 0.9*v_low, 1.1*v_high,
           dir_fig+"Line_ratio_N2Ha", cmap='rainbow', **pltFlags)

# Error map
plot_2Dmap(plt_Data/e_plt_Data, r"${\rm [NII]\lambda 6584/H\alpha}$ flux ratio S/N map", 0.0, 3.0,
           dir_fig+"Line_ratio_eN2Ha", cmap='rainbow', **pltFlags)
# ----- END: [NII]6584 / H alpha flux ratio map ----- #


# ----- START: Oxygen abundance map (N2 method) ----- #
N2 = np.log10(plt_Data)
logOH = 8.743 + 0.462*N2
v_low, v_high = np.percentile(logOH[np.isnan(logOH) == False], [1.0, 99.0])
plot_2Dmap(logOH, r"${\rm 12+log(O/H)}$ map (N2 method)", 8.1, 8.7,
           dir_fig+"Map_logOH_N2", cmap='rainbow', **pltFlags)

def cal_fwm_logOH(regfile, metal_data):
    metal = copy.deepcopy(metal_data)
    metal[metal <= 0.] = np.nan
    msk = (np.isnan(metal) == True)

    x0, y0, a, b, pa = read_region(regfile, regtype='ellipse')
    x0, y0 = x0[0]-1, y0[0]-1
    a, b = a[0], b[0]
    pa = pa[0]*np.pi/180.

    ap = EAp(positions=[(x0, y0)], a=a, b=b, theta=pa)

    tbl0 = apphot(data=flx_Data, apertures=ap, mask=msk)
    flx_sum_disk = tbl0['aperture_sum'].data[0]
    flx_sum_tail = np.sum(flx_Data[~msk]) - flx_sum_disk

    tbl1 = apphot(data=flx_Data*metal, apertures=ap, mask=msk)
    logOH_fwm_disk = tbl1['aperture_sum'].data[0] / flx_sum_disk
    logOH_fwm_tail = (np.sum(flx_Data[~msk]*metal[~msk]) - tbl1['aperture_sum'].data[0]) / flx_sum_tail

    return [logOH_fwm_disk, logOH_fwm_tail]

# Flux-weighted mean
val_N2 = (np.isnan(N2) == False)
fwm_logOH = weighted_mean(data=logOH[val_N2], weights=Halpha_flux_2D[val_N2])
print(f"Flux-weighted mean of log O/H (N2 method) : {fwm_logOH:.3f}")

# # 1. HST
# logOH_disk_hst, logOH_tail_hst = cal_fwm_logOH("HST_boundary_1sig_transformed.reg", logOH)
# print(f"log O/H (HST, N2 method) : disk - {logOH_disk_hst:.3f}, tail - {logOH_tail_hst:.3f}")

# 2. Gemini
logOH_disk_gem, logOH_tail_gem = cal_fwm_logOH("GMOS_boundary_1sig.reg", logOH)
print(f"log O/H (Gemini, N2 method) : disk - {logOH_disk_gem:.3f}, tail - {logOH_tail_gem:.3f}")
# ----- END: Oxygen abundance map (N2 method) ----- #


# ----- START: O3N2 map ----- #
O3N2 = np.log10((OIII5007_flux_2D / Hbeta_flux_2D) * (Halpha_flux_2D / NII6584_flux_2D))
zero2_cnd = (zero_cnd | (OIII5007_snr_2D < 3.0) | (Hbeta_snr_2D < 3.0) | (NII6584_snr_2D < 3.0))
O3N2[zero2_cnd] = 0.
O3N2[((O3N2 == 0.) | (np.isinf(O3N2) == True))] = np.nan
v_low, v_high = np.percentile(O3N2[np.isnan(O3N2) == False], [1.0, 99.0])
plot_2Dmap(O3N2, "O3N2 map", v_low, v_high, dir_fig+"Map_O3N2", cmap='rainbow', **pltFlags)
# ----- END: O3N2 map ----- #


# ----- START: Oxygen abundance map (O3N2) ----- #
logOH2 = 8.533-0.214*O3N2
v_low, v_high = np.percentile(logOH2[np.isnan(logOH2) == False], [1.0, 99.0])
plot_2Dmap(logOH2, r"${\rm 12+log(O/H)}$ map (O3N2 method)", 8.1, 8.7,
           dir_fig+"Map_logOH_O3N2", cmap='rainbow', **pltFlags)

# Flux-weighted mean
val_O3N2 = (np.isnan(O3N2) == False)
fwm_logOH2 = weighted_mean(data=logOH2[val_O3N2], weights=Halpha_flux_2D[val_O3N2])
print(f"Flux-weighted mean of log O/H (O3N2 method) : {fwm_logOH2:.3f}")

# # 1. HST
# logOH2_disk_hst, logOH2_tail_hst = cal_fwm_logOH("HST_boundary_1sig_transformed.reg", logOH2)
# print(f"log O/H (HST, O3N2 method) : disk - {logOH2_disk_hst:.3f}, tail - {logOH2_tail_hst:.3f}")

# 2. Gemini
logOH2_disk_gem, logOH2_tail_gem = cal_fwm_logOH("GMOS_boundary_1sig.reg", logOH2)
print(f"log O/H (Gemini, O3N2 method) : disk - {logOH2_disk_gem:.3f}, tail - {logOH2_tail_gem:.3f}")
# ----- END: Oxygen abundance map (O3N2) ----- #


# ----- START: [OIII]5007 / H beta flux ratio map ----- #
plt_Data, e_plt_Data = get_line_ratio(OIII5007_flux_2D, Hbeta_flux_2D, OIII5007_snr_2D, Hbeta_snr_2D,
                                      e_OIII5007_flux_2D, e_Hbeta_flux_2D)
v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
plot_2Dmap(plt_Data, r"${\rm [OIII]\lambda 5007/H\beta}$ flux ratio map", 0.9*v_low, 1.1*v_high,
           dir_fig+"Line_ratio_O3Hb", cmap='rainbow', **pltFlags)
# ----- END: [OIII]5007 / H beta flux ratio map ----- #


# ----- START: [SII] / H alpha flux ratio map ----- #
S2_flux = SII6717_flux_2D + SII6731_flux_2D
S2_snr = S2_flux / np.sqrt(e_SII6717_flux_2D**2 + e_SII6731_flux_2D**2)
plt_Data, e_plt_Data = get_line_ratio(S2_flux, Halpha_flux_2D, S2_snr, Halpha_snr_2D,
                                      np.sqrt(e_SII6717_flux_2D**2 + e_SII6731_flux_2D**2),
                                      e_Halpha_flux_2D)
v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
plot_2Dmap(plt_Data, r"${\rm [SII]\lambda\lambda 6717,6731/H\alpha}$ flux ratio map",
           0.9*v_low, 1.1*v_high, dir_fig+"Line_ratio_S2Ha", cmap='rainbow', **pltFlags)
# ----- END: [SII]6731 / H alpha flux ratio map ----- #


# ----- START: [SII]6717 / [SII]6731 flux ratio map ----- #
plt_Data, e_plt_Data = get_line_ratio(SII6717_flux_2D, SII6731_flux_2D, SII6717_snr_2D, SII6731_snr_2D,
                                      e_SII6717_flux_2D, e_SII6731_flux_2D)
v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
plot_2Dmap(plt_Data, r"${\rm [SII]\lambda 6717/[SII]\lambda 6731}$ flux ratio map", 0.9*v_low, 1.1*v_high,
           dir_fig+"Line_ratio_S2S2", cmap='rainbow', **pltFlags)
# ----- END: [SII]6717 / [SII]6731 flux ratio map ----- #


# ----- START: [NII]6584 / [NII]6548 flux ratio map ----- #
plt_Data, e_plt_Data = get_line_ratio(NII6584_flux_2D, NII6548_flux_2D, NII6584_snr_2D, NII6548_snr_2D,
                                      e_NII6584_flux_2D, e_NII6548_flux_2D)
v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
plot_2Dmap(plt_Data, r"${\rm [NII]\lambda 6584/[NII]\lambda 6548}$ flux ratio map", 0.9*v_low, 1.1*v_high,
           dir_fig+"Line_ratio_N2N2", cmap='rainbow', **pltFlags)
# ----- END: [NII]6584 / [NII]6548 flux ratio map ----- #


# ----- Saving the results ----- #
np.savez('plot_data.npz', sflx=sflx, rvd=Rvd, vdd=Vdd, sfrd=SFRD)
df_res1 = pd.Series(data = {"FW-mean vdisp": fwm_vdisp,
                            "FW-mean Hab": fwm_Hab,
                            "FW-mean Hab err": e_fwm_Hab,
                            "SFR total": np.sum(SFR[val_SFR]),
                            "SFR total err": np.sqrt(np.sum(e_SFR[val_SFR]**2)),
                            "Mean A_V": np.mean(A_V[val_SFR]),
                            "SFR disk": SFR_disk_gem,
                            "SFR disk err": e_SFR_disk_gem})


### Additional analysis

# ----- START: BPT diagram ----- #
x_Dat = np.log10(NII6584_flux_2D/Halpha_flux_2D)
y_Dat = np.log10(OIII5007_flux_2D/Hbeta_flux_2D)

snr2_cnd = ((NII6584_snr_2D < 3.0) | (Halpha_snr_2D < 3.0) | \
           (OIII5007_snr_2D < 3.0) | (Hbeta_snr_2D < 3.0))
zero2_cnd = (zero_cnd | snr2_cnd)

x_Dat[zero2_cnd] = np.nan
y_Dat[zero2_cnd] = np.nan

BPT_SFG = (y_Dat < 1.3+0.61/x_Dat)
BPT_comp = ((y_Dat > 1.3+0.61/x_Dat) & (y_Dat < 1.19+0.61/(x_Dat-0.47)))
BPT_AGN = ((y_Dat > 1.19+0.61/(x_Dat-0.47)) & \
           (y_Dat-1.4899 > ((0.2949- -0.2135)/(1.4899-0.4548))*(x_Dat-0.2949)))
BPT_LINER = ((y_Dat > 1.19+0.61/(x_Dat-0.47)) & \
             (y_Dat-1.4899 < ((0.2949- -0.2135)/(1.4899-0.4548))*(x_Dat-0.2949)))

x_Dat_flat = x_Dat.flatten()
y_Dat_flat = y_Dat.flatten()

fig, ax = plt.subplots(1, 1, figsize=(8,5))
plt.suptitle("BPT diagram", x=0.5, ha='center', y=0.96, va='top', fontsize=20.0)
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
plt.suptitle("The spatial map of the BPT class", x=0.5, ha='center', y=0.96, va='top', fontsize=20.0)
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
            rect = plt.Rectangle((-pixel_scale*bpt_data.shape[1]/2 + pixel_scale*np.argwhere(bpt_data == True)[j,1],
                                  pixel_scale*bpt_data.shape[0]/2 - pixel_scale*np.argwhere(bpt_data == True)[j,0]),
                                 pixel_scale, pixel_scale, facecolor=cmap(c_id)[:-1], alpha=0.8)
            ax.add_patch(rect)

ax.contour(X_coord, Y_coord[::-1], sflx, levels=lvs, linewidths=lws, colors=cs, alpha=0.6)
p0, = ax.plot(-100.0, -100.0, '-', linewidth=2.5, color='gray', alpha=0.6,
              label=r"H${\rm \alpha}$ flux contour")
ax.legend(handles=[p0], fontsize=13.0, loc='lower left',
          handlelength=2.5, frameon=True, borderpad=0.8,
          framealpha=0.8, edgecolor='gray')

# The orientations
ax.arrow(pltFlags['x0']+pltFlags['sign']*0.025, pltFlags['y0'],
         pltFlags['L']*np.sin(pltFlags['theta0']), pltFlags['L']*np.cos(pltFlags['theta0']),
         width=0.06, head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.arrow(pltFlags['x0'], pltFlags['y0']+pltFlags['sign']*0.025,
         -pltFlags['L']*np.cos(pltFlags['theta0']), pltFlags['L']*np.sin(pltFlags['theta0']),
         width=0.06, head_width=0.18, head_length=0.18, fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(pltFlags['xE'], pltFlags['yE'], 'E', fontsize=15.0, fontweight='bold', color='blueviolet')
ax.text(pltFlags['xN'], pltFlags['yN'], 'N', fontsize=15.0, fontweight='bold', color='blueviolet')

# Scale bar
kpc5 = 5.0 / ang_scale
ax.arrow(2.0, -1.85, kpc5, 0., width=0.07, head_width=0., head_length=0.,
          fc='blueviolet', ec='blueviolet', alpha=0.9)
ax.text(2.1, -2.2, '5 kpc', fontsize=15.0, fontweight='bold', color='blueviolet')

plt.savefig(dir_fig+'BPT_map.pdf')
plt.savefig(dir_fig+'BPT_map.png', dpi=300)
plt.close()
# ----- END: BPT spatial map ----- #

# '''
# # ----- Applying NebulaBayes ----- #
# import NebulaBayes 
# dir_NB = "/".join(os.path.abspath(NebulaBayes.__file__).split("/")[:-1])
# from NebulaBayes import NB_Model

# # These HII-region optical emission-line fluxes have already been dereddened
# norm_line = "Halpha"
# linelist_NBinput = ["Hbeta", "OII3729", "OIII5007", "Halpha", "NII6583", "SII6716", "SII6731"]
# linelist_GEMname = ["Hbeta", "OII3727", "OIII5007", "Halpha", "NII6584", "SII6717", "SII6731"]

# df_ll = pd.read_csv(dir_NB+"/grids/Linelist.csv")
# wavlist_NBinput = []
# for l in linelist_NBinput:
#     wavlist_NBinput.append(df_ll['Lambda_AA'][df_ll['Grid_name'] == l].values[0])

# # Set outputs:
# OUT_DIR = dir_fig+"NB_HII"
# if (glob.glob(OUT_DIR) == []):
#     os.system("mkdir "+OUT_DIR)
# else:
#     os.system("rm -rf "+OUT_DIR+"/*")

# # Initialize the NB_Model, which loads and interpolates the model flux grids:
# grid_table_file = os.path.join(dir_NB, "grids", "NB_HII_grid.fits.gz")
# BinTableHDU_0 = fits.getdata(grid_table_file, ext=0)
# DF_grid = Table(BinTableHDU_0).to_pandas()
# DF_grid = DF_grid[DF_grid["log P/k"] == 5.8]    # Fixing log P/k = 5.8
# grid_params = ["log U", "12 + log O/H"]

# Ngrid = (100, 200)    # less than 60000
# NB_Model_HII = NB_Model(DF_grid, grid_params, line_list=linelist_NBinput,
#                         interpd_grid_shape=Ngrid, grid_error=0.1)

# # Setting custom prior
# Result0_HII = NB_Model_HII([1.]*len(linelist_NBinput), [0.1]*len(linelist_NBinput),
#                            linelist_NBinput, norm_line=norm_line)
# grid_spec = Result0_HII.Posterior.Grid_spec
# param_cut = [[-3.75, -3.0], [7.663, 8.76]]
# prior = np.ones(grid_spec.shape)
# prior_func = "Gaussian"    # "Uniform", "Gaussian"
# sigma_func = [0.1, 0.1]
# for p in np.arange(len(grid_params)):
#     param_idx = grid_spec.param_names.index(grid_params[p])
#     all_values = grid_spec.param_values_arrs[param_idx]    # Sorted 1D array

#     cut_lo, cut_hi = param_cut[p][0], param_cut[p][1]
#     prior_1D = np.ones_like(all_values)

#     if (prior_func == "Uniform"):
#         prior_1D[all_values >= cut_hi] = 0.
#         prior_1D[all_values <= cut_lo] = 0.

#     if (prior_func == "Gaussian"):
#         gauss_lo = np.exp(-((all_values-cut_lo)/sigma_func[p])**2 / 2.)
#         gauss_hi = np.exp(-((all_values-cut_hi)/sigma_func[p])**2 / 2.)
#         prior_1D[all_values <= cut_lo] = gauss_lo[all_values <= cut_lo]
#         prior_1D[all_values >= cut_hi] = gauss_hi[all_values >= cut_hi]

#     # This array has only one dimension.  Construct a slice to use numpy
#     # "broadcasting" to apply the 1D prior over the whole 2D grid:
#     # (The following method works in nD, although here it's a 2D grid)
#     slice_NLR = [np.newaxis for _ in grid_spec.shape]
#     slice_NLR[param_idx] = slice(None)  # "[slice(None)]" means "[:]"
#     slice_NLR = tuple(slice_NLR)
#     prior *= prior_1D[slice_NLR]

# # Initializing parameters
# NB_logOH_2D = np.zeros_like(Halpha_flux_2D)
# NB_logU_2D = np.zeros_like(Halpha_flux_2D)
# fault_bin = []

# # Calculating weight-mean flux ratio in advance
# flx0_Data = copy.deepcopy(Halpha_flux_2D)
# e_flx0_Data = copy.deepcopy(e_Halpha_flux_2D)
# flx0_sum = np.sum(flx0_Data)
# e_flx0_sum = np.sqrt(np.sum(e_flx0_Data**2.))

# wm_ratios, e_wm_ratios = [], []
# for l in linelist_GEMname:
#     exec("flx1_Data = copy.deepcopy("+l+"_flux_2D)")
#     exec("e_flx1_Data = copy.deepcopy(e_"+l+"_flux_2D)")
#     exec("snr_Data = copy.deepcopy("+l+"_snr_2D)")
    
#     flx_ratio = flx1_Data / flx0_Data
#     snr2_cnd = (snr_Data < 3.0)
#     zero2_cnd = (zero_cnd | snr2_cnd)
#     flx_ratio[zero2_cnd] = 0.
#     flx_ratio[flx_ratio == 0.] = np.nan
#     flx_ratio[np.isinf(flx_ratio) == True] = np.nan
#     e_flx_ratio = flx_ratio * np.sqrt((e_flx0_Data/flx0_Data)**2 + (e_flx1_Data/flx1_Data)**2)
#     val = (np.isnan(flx_ratio) == False)
#     wm, e_wm = weighted_mean(flx_ratio[val], flx0_Data[val], e_data=e_flx_ratio[val], e_weights=e_flx0_Data[val])
#     wm_ratios.append(wm)
#     e_wm_ratios.append(e_wm)

# # Running NebulaBayse for each Voronoi bin
# snr2_cnd = (NII6584_snr_2D < 3.0)
# zero2_cnd = (zero_cnd | snr2_cnd)

# for i in np.arange(nvbin):
#     if (np.unique(zero2_cnd[data_vbin == i])[0] == False):
#         print(f"\n----- Running NebulaBayes for bin {i:d} -----\n")
#         obs_fluxes, obs_errs = [], []
#         for l in linelist_GEMname:
#             idx_l = linelist_GEMname.index(l)
#             exec("fl = "+l+"_flux_2D[data_vbin == i]")
#             exec("e_fl = e_"+l+"_flux_2D[data_vbin == i]")
#             fl_input, e_fl_input = np.unique(fl)[0], np.unique(e_fl)[0]
#             if (fl_input <= 0.0):
#                 fl0 = np.unique(Halpha_flux_2D[data_vbin == i])[0]
#                 e_fl0 = np.unique(e_Halpha_flux_2D[data_vbin == i])[0]
#                 fl_input = fl0 * wm_ratios[idx_l]
#                 e_fl_input = fl_input * np.sqrt((e_fl0/fl0)**2 + (e_wm_ratios[idx_l]/wm_ratios[idx_l])**2)
#             obs_fluxes.append(fl_input)
#             obs_errs.append(e_fl_input)

#         kwargs = {"norm_line": norm_line, "deredden": False,
#                   "obs_wavelengths": wavlist_NBinput, "prior": prior,
#                   # "prior": "Uniform",
#                   "prior_plot": os.path.join(OUT_DIR, f"1_HII_prior_plot_bin{i:d}.pdf"),
#                   "likelihood_plot": os.path.join(OUT_DIR, f"1_HII_likelihood_plot_bin{i:d}.pdf"),
#                   "posterior_plot": os.path.join(OUT_DIR, f"1_HII_posterior_plot_bin{i:d}.pdf"),
#                   "estimate_table": os.path.join(OUT_DIR, f"1_HII_param_estimates_bin{i:d}.csv"),
#                   "best_model_table": os.path.join(OUT_DIR, f"1_HII_best_model_bin{i:d}.csv")
#                   }

#         try:
#             Result_HII = NB_Model_HII(obs_fluxes, obs_errs, linelist_NBinput, **kwargs)
#             Estimate_table = Result_HII.Posterior.DF_estimates  # pandas DataFrame
#             print("\nParameter estimate table:")
#             print(Estimate_table)

#             for p in grid_params:
#                 est = Estimate_table.loc[p, "Estimate"]
#                 low = Estimate_table.loc[p, "CI68_low"]
#                 high = Estimate_table.loc[p, "CI68_high"]
#                 print("\nThe measured "+p+" = "
#                       "{0:.2f}^{{+{1:.2f}}}_{{-{2:.2f}}}".format(est, *(high-est, est-low)))

#             best_model_dict = Result_HII.Posterior.best_model
#             print("\nBest model table:")
#             print(best_model_dict["table"])  # pandas DataFrame

#             NB_logOH_2D[data_vbin == i] = Estimate_table.loc["12 + log O/H", "Estimate"]
#             NB_logU_2D[data_vbin == i] = Estimate_table.loc["log U", "Estimate"]
        
#         except ValueError:
#             # print(ValueError)
#             continue
#     else:
#         print(f"\n----- NebulaBayes is not available for bin {i:d} -----\n")
#         fault_bin.append(i)


# # ----- START: Oxygen abundance map (NebulaBayes, HII) ----- #
# plt_Data = copy.deepcopy(NB_logOH_2D)
# plt_Data[plt_Data == 0] = np.nan
# v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
# plot_2Dmap(plt_Data, r"${\rm 12+log(O/H)}$ map (NebulaBayes, HII)", 8.1, 8.7,
#            dir_fig+"Map_logOH_NBH2", cmap='rainbow')

# # Flux-weighted mean
# val_NB = (np.isnan(plt_Data) == False)
# fwm_logOH3 = weighted_mean(data=plt_Data[val_NB], weights=Halpha_flux_2D[val_NB])
# print(f"Flux-weighted mean of log O/H (NB, HII) : {fwm_logOH3:.3f}")

# # # 1. HST
# # logOH2_disk_hst, logOH2_tail_hst = cal_fwm_logOH("HST_boundary_1sig_transformed.reg", logOH2)
# # print(f"log O/H (HST, O3N2 method) : disk - {logOH2_disk_hst:.3f}, tail - {logOH2_tail_hst:.3f}")

# # 2. Gemini
# logOH3_disk_gem, logOH3_tail_gem = cal_fwm_logOH("GMOS_boundary_1sig.reg", plt_Data)
# print(f"log O/H (Gemini, NB, HII) : disk - {logOH3_disk_gem:.3f}, tail - {logOH3_tail_gem:.3f}")
# # ----- END: Oxygen abundance map (NebulaBayes, HII) ----- #


# # ----- START: Ionization parameter map ----- #
# plt_Data = copy.deepcopy(NB_logU_2D)
# plt_Data[plt_Data == 0] = np.nan
# plt_Data += np.log10(c*1.0e+5)
# v_low, v_high = np.percentile(plt_Data[np.isnan(plt_Data) == False], [1.0, 99.0])
# plot_2Dmap(plt_Data, r"${\rm log(q)}$ map (NebulaBayes, HII)", 6.8, 7.5,
#            dir_fig+"Map_logQ_NBH2", cmap='rainbow')
# # ----- END: Ionization parameter map ----- #


# # ----- Saving the results ----- #
# np.savez('plot_data2.npz', logOH_1=NB_logOH_2D, logU_1=NB_logU_2D)
# df_res2 = pd.Series(data = {"N2_total": fwm_logOH,
#                             "N2_disk": logOH_disk_gem,
#                             "N2_tail": logOH_tail_gem,
#                             "O3N2_total": fwm_logOH2,
#                             "O3N2_disk": logOH2_disk_gem,
#                             "O3N2_tail": logOH2_tail_gem,
#                             "NBH2_total": fwm_logOH3,
#                             "NBH2_disk": logOH3_disk_gem,
#                             "NBH2_tail": logOH3_tail_gem})
# '''

# Printing the running time
print('--- %.4f seconds ---' %(time.time()-start_time))
