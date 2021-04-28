#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 16:29:33 2019

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import glob, os, copy
import g0_init_cfg as ic
from astropy.io import fits
from reg_saoimage import read_region


# ----- Importing IRAF from the root directory ----- #
current_dir = os.getcwd()
os.chdir(ic.dir_iraf)

from pyraf import iraf
from pyraf.iraf import gemini, gmos

os.chdir(current_dir)
iraf.chdir(current_dir)


# ---------- Fixing a bad column ---------- #
std = np.loadtxt(ic.lst_std, dtype=str)
if (std.size > 1):
    raise ValueError("Please check if there is only one image for the standard star.")
std0 = std.item(0)

# 2-slit mode
dt1, hd1 = fits.getdata('eqxbrg'+std0+'.fits', ext=2, header=True)
if (ic.nslit == 2):
	dt2, hd2 = fits.getdata('eqxbrg'+std0+'.fits', ext=5, header=True)


# ----- Reading masking coordinates ----- #
if (ic.nslit == 1):
    fix_std = {'slit1': [[], []], 'slit2': [[], []]}
if (ic.nslit == 2):
    fix_std = {'slit1': [[], []], 'slit2': [[], []]}

for i in 1+np.arange(ic.nslit):
    regfile = "badcol/std1_slit{0:1d}.reg".format(i)
    if (glob.glob(regfile) == []):
        continue
    else:
        x0, y0, xs, ys, _ = read_region(regfile, regtype="box")
        fix_std["slit{0:1d}".format(i)][0].append(int(round(x0[0]))-1)
        fix_std["slit{0:1d}".format(i)][1].append(int(round(xs[0]))-1)


for i in np.arange(len(iraf.type(ic.lst_std, Stdout=1))):
    std = iraf.type(ic.lst_std, Stdout=1)[i].strip()
    exec('dic = copy.deepcopy(fix_sci{0:02d})'.format(i+1))
    for j in np.arange(ic.nslit):
    	fmsk = 'mskbadcol_'+std0+'_{0:1d}.txt'.format(j+1)
    	s = dic['slit{0:1d}'.format(j+1)]
    	exec('d = dt{0:1d}'.format(j+1))
    	for k in np.arange(len(s[0])):
    		if (k == 0):
    			com = "echo '{0:d} {1:d} 1 {2:d}' > ".format(s[0][k]-s[1][k]/2, s[0][k]+s[1][k]/2, d.shape[0])+fmsk
    			os.system(com)
    		else:
    			com = "echo '{0:d} {1:d} 1 {2:d}' >> ".format(s[0][k]-s[1][k]/2, s[0][k]+s[1][k]/2, d.shape[0])+fmsk
    			os.system(com)
    	iraf.text2mask(fmsk, fmsk.strip('.txt')+'.pl', d.shape[1], d.shape[0])

iraf.imdelete('xeqxbrg@'+ic.lst_std)
iraf.imdelete('tmp@'+ic.lst_std)
iraf.imdelete('tmpdq@'+ic.lst_std)

for std in iraf.type(ic.lst_std, Stdout=1):
    std = std.strip()
    iraf.copy('eqxbrg'+std0+'.fits', 'tmp'+std0+'.fits')
    for j in np.arange(ic.nslit):
    	iraf.proto.fixpix('tmp'+std0+'.fits[sci,{0:1d}]'.format(j+1),
    		              'mskbadcol_'+std0+'_{0:1d}.pl'.format(j+1),
    		              linterp='1,2,3,4')
    iraf.copy('tmp'+std0+'.fits', 'xeqxbrg'+std0+'.fits')
    for j in np.arange(ic.nslit):
    	iraf.imarith('mskbadcol_'+std0+'_{0:1d}.pl'.format(j+1), '+',
    		         'xeqxbrg'+std0+'.fits[dq,{0:1d}]'.format(j+1),
    		         'tmpdq'+std0+'.fits[dq,{0:1d}]'.format(j+1))
    	iraf.imcopy('tmpdq'+std0+'.fits[dq,{0:1d}][*,*]'.format(j+1),
    		        'xeqxbrg'+std0+'.fits[dq,{0:1d}][*,*]'.format(j+1))
    	iraf.proto.fixpix('tmp'+std0+'.fits[sci,{0:1d}]'.format(j+1),
    		              'mskbadcol_'+std0+'_{0:1d}.pl'.format(j+1),
    		              linterp='1,2,3,4')    


for std in iraf.type(ic.lst_std, Stdout=1):
	for j in np.arange(ic.nslit):
		os.system('ds9 &')
		iraf.sleep(5.0)
		iraf.display(image = 'eqxbrg'+std0+'.fits[sci,{0:1d}]'.format(j+1), frame=1)
		iraf.display(image = 'xeqxbrg'+std0+'.fits[sci,{0:1d}]'.format(j+1), frame = 2)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
