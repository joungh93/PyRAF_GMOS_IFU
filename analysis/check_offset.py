#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 15:44:25 2020

@author: jlee
"""

import time
start_time = time.time()

import numpy as np
import glob, os
import init_cfg as ic
from astropy.io import fits
from scipy import ndimage


# ----- Loading Ha sum images ----- #
current_dir = os.getcwd()
working_dir = ic.dir_cmb
os.chdir(working_dir)
# Ha_list = sorted(glob.glob("Ha_sum-*.fits"))


# ----- Calculating WCS offsets ----- #

# WCS reference
idx_ref = np.where([ic.cube_ref in cb for cb in ic.cube_list])[0][0]
h01 = fits.getheader(ic.cube_list[idx_ref], ext=0)
PA = h01['PA']*np.pi/180.

# Writing WCS offset file
f = open('offset.txt', 'w')
for i in np.arange(len(ic.cube_list)):
	h = fits.getheader(ic.cube_list[i], ext=0)

	offset_RA = (h['CRVAL1']-h01['CRVAL1'])*np.cos(h01['CRVAL2'])*3600.0    # arcsec
	offset_Dec = (h['CRVAL2']-h01['CRVAL2'])*3600.0    # arcsec

	offset_X = (+offset_RA*np.cos(PA) - offset_Dec*np.sin(PA)) / ic.pixel_scale    # pixel
	offset_Y = (+offset_RA*np.sin(PA) + offset_Dec*np.cos(PA)) / ic.pixel_scale    # pixel

	print(ic.cube_name[i]+f" - RA offset: {offset_RA:.3f}, Dec offset: {offset_Dec:.3f}")
	print(f"                 X offset: {offset_X:.3f} pix, Y offset: {offset_Y:.3f}")
	print(f"                 Reference pixel ({h['CRPIX1']:.3f}, {h['CRPIX2']:.3f}) \n")
	f.write(f"{offset_X:.3f}  {offset_Y:.3f} \n")
f.close()


# ----- Shifting & combining Ha sum images ----- #

# Running an initial shift task
offset_X, offset_Y = np.loadtxt('offset.txt').T
dhs0 = fits.getdata(sorted(glob.glob("Ha_sum-*.fits"))[0], ext=0, header=False)
dhs2 = np.zeros((len(ic.cube_list), dhs0.shape[0], dhs0.shape[1]))
if (len(ic.cube_spa_off) > 0):
	dhs3 = np.zeros((len(ic.cube_list)-len(ic.cube_spa_off), dhs0.shape[0], dhs0.shape[1]))
j = 0
for i in np.arange(len(ic.cube_list)):
	dhs, hdr = fits.getdata("Ha_sum-"+ic.cube_name[i]+".fits", ext=0, header=True)
	dhs_shifted = ndimage.shift(dhs, shift=(-offset_Y[i], -offset_X[i]), mode='nearest')
	fits.writeto('al1_Ha_sum-'+ic.cube_name[i]+'.fits', dhs_shifted, hdr, overwrite=True)
	dhs2[i,:,:] = dhs_shifted
	if ((len(ic.cube_spa_off) > 0) & (ic.cube_name[i] not in ic.cube_spa_off)):
		dhs3[j,:,:] = dhs_shifted
		j += 1

# Combining the shifted images
cdhs = np.median(dhs2, axis=0)

##### Needed to be revised manually after visual check! #####
if ((len(ic.cube_spa_off) > 0) & (ic.overwrite == True)):
	cdhs[:,61:] = np.median(dhs3[:,:,61:], axis=0)
##########

fits.writeto('al1_fcomb.fits', cdhs, overwrite=True)
disp = "ds9 -scalemode zscale -scale lock yes -frame lock image "
tile = "-tile grid mode manual -tile grid layout "
cross = "-mode crosshair -lock crosshair image "
os.system(disp+cross+"Ha_sum-*.fits al1_Ha_sum-*.fits "+tile+f"{len(ic.cube_list):d} 2 &")
os.system(disp+"al1_fcomb.fits &")


# # ----- Running IRAF/xregister task for the shifted Ha sum images ----- #
# iraf.images()
# iraf.immatch()
# iraf.unlearn('xregister')

# os.system('ls -1 al1_Ha_sum-*.fits > input_al1.lis')
# os.system('rm -rfv shifts.db rgal1_Ha_sum-*.fits')
# iraf.xregister(input='@input_al1.lis', reference='al1_Ha_sum-N20190613S0237.fits', regions='[42:62, 10:36]',
#                shifts='shifts.db', output='rg//@input_al1.lis', databasefmt='no', correlation='fourier',
#                interp_type='linear', boundary_type='nearest')

# rgal1_list = sorted(glob.glob('rgal1_Ha_sum*.fits'))

# f1 = open('rgal1_comb1.lis', 'w')
# f2 = open('rgal1_comb2.lis', 'w')
# for i in np.arange(len(ic.cube_list)):
# 	cb = ic.cube_list[i].split('/')[-1].split('cstxeqxbrg')[-1].split('_3D.fits')[0]
# 	# cb = ic.cube_list[i].split('_3D.fits')[0].split(ic.dir_redux+'cstxeqxbrg')[1]
# 	if (cb not in ic.cube_spa_off):
# 		f1.write(rgal1_list[i]+'\n')
# 	else:
# 		f2.write(rgal1_list[i]+'\n')
# f1.close()
# f2.close()

# # Combining w/ IRAF/imcombine task
# outimg = ['rgal1_Ha_comb1.fits', 'rgal1_Ha_comb2.fits']
# for i in np.arange(len(outimg)):
# 	os.system('rm -rfv '+outimg[i])
# 	iraf.imcombine(input='@rgal1_comb{0:1d}.lis'.format(i+1), output=outimg[i], combine='median')
# 	exec('dat_rgal1_comb{0:1d} = fits.getdata(outimg[i])'.format(i+1))

# dat_rgal1_fcomb = 0.5*(dat_rgal1_comb1+dat_rgal1_comb2)
# dat_rgal1_fcomb[:, 62-1:] = dat_rgal1_comb1[:,62-1:]

# fits.writeto('rgal1_fcomb.fits', dat_rgal1_fcomb, overwrite=True)

# os.system('ds9 -scalemode zscale -scale lock yes -frame lock image rgal1_fcomb.fits &')


os.chdir(current_dir)


# Printing the running time
print('--- %.4f seconds ---' %(time.time()-start_time))
