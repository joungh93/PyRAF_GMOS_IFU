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


# ----- Loading Ha sum images ----- #
current_dir = os.getcwd()
working_dir = ic.dir_cmb
os.chdir(working_dir)
Ha_list = sorted(glob.glob("Ha_sum-*.fits"))


# ----- Calculating WCS offsets ----- #

# WCS reference
h01 = fits.getheader(ic.cube_list[0], ext=0)
PA = h01['PA']*np.pi/180.

# Writing WCS offset file
f = open('offset.txt', 'w')
for i in np.arange(len(ic.cube_list)):
	h = fits.getheader(ic.cube_list[i], ext=0)

	offset_RA = (h['CRVAL1']-h01['CRVAL1'])*np.cos(h01['CRVAL2'])*3600.0    # arcsec
	offset_Dec = (h['CRVAL2']-h01['CRVAL2'])*3600.0    # arcsec

	offset_X = (-offset_RA*np.cos(PA) + offset_Dec*np.sin(PA)) / ic.pixel_scale    # pixel
	offset_Y = (offset_RA*np.sin(PA) + offset_Dec*np.cos(PA)) / ic.pixel_scale    # pixel

	print(ic.cube_name[i]+f" - RA offset: {offset_RA:.3f}, Dec offset: {offset_Dec:.3f}")
	print(f"                 X offset: {offset_X:.3f} pix, Y offset: {offset_Y:.3f}")
	f.write(f"{offset_X:.3f}  {offset_Y:.3f} \n")
f.close()


# # ----- Shifting & combining Ha sum images ----- #

# # Shifting w/ IRAF/imshift task
# iraf.images()
# iraf.imgeom()
# iraf.unlearn('imshift')

# os.system('ls -1 Ha_sum-*.fits > input_Ha_sum.lis')
# os.system('rm -rfv al1_Ha_sum*.fits')
# iraf.imshift(input='@input_Ha_sum.lis', output='al1_//@input_Ha_sum.lis', shifts_file='offset.txt',
#              interp_type='linear', boundary_type='nearest')

# al1_list = sorted(glob.glob('al1_Ha_sum*.fits'))

# f1 = open('al1_comb1.lis', 'w')
# f2 = open('al1_comb2.lis', 'w')
# for i in np.arange(len(ic.cube_list)):
# 	cb = ic.cube_list[i].split('/')[-1].split('cstxeqxbrg')[-1].split('_3D.fits')[0]
# 	# cb = ic.cube_list[i].split('_3D.fits')[0].split(ic.dir_redux+'cstxeqxbrg')[1]
# 	if (cb not in ic.cube_spa_off):
# 		f1.write(al1_list[i]+'\n')
# 	else:
# 		f2.write(al1_list[i]+'\n')
# f1.close()
# f2.close()

# # Combining w/ IRAF/imcombine task (for reference)
# outimg = ['al1_Ha_comb1.fits', 'al1_Ha_comb2.fits']
# for i in np.arange(len(outimg)):
# 	os.system('rm -rfv '+outimg[i])
# 	iraf.imcombine(input='@al1_comb{0:1d}.lis'.format(i+1), output=outimg[i], combine='median')
# 	exec('dat_al1_comb{0:1d} = fits.getdata(outimg[i])'.format(i+1))

# dat_al1_fcomb = 0.5*(dat_al1_comb1+dat_al1_comb2)
# dat_al1_fcomb[:, 62-1:] = dat_al1_comb1[:,62-1:]

# fits.writeto('al1_fcomb.fits', dat_al1_fcomb, overwrite=True)


# os.system('ds9 -scalemode zscale -scale lock yes -frame lock image al1_fcomb.fits &')


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
print('--- %s seconds ---' %(time.time()-start_time))
