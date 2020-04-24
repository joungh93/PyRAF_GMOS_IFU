#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 16:29:33 2019

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import glob, os
import g0_init_cfg as ic


# ----- Importing IRAF from the root directory ----- #
current_dir = os.getcwd()
os.chdir(ic.dir_iraf)

from pyraf import iraf
from pyraf.iraf import gemini, gmos

os.chdir(current_dir)
iraf.chdir(current_dir)


# ---------- Fixing a bad column ---------- #
from astropy.io import fits

for sci in iraf.type(ic.lst_sci, Stdout=1):
    sci = sci.strip()
    fits.open('eqxbrg'+sci+'.fits').info()
    # 2-slit mode
    dt1, hd1 = fits.getdata('eqxbrg'+sci+'.fits', ext=2, header=True)
    dt2, hd2 = fits.getdata('eqxbrg'+sci+'.fits', ext=5, header=True)

    # os.system('ds9 &')
    # iraf.sleep(5.0)
    # iraf.display(image = 'eqxbrg'+sci+'.fits[sci,1]', frame = 1)
    # iraf.display(image = 'eqxbrg'+sci+'.fits[sci,2]', frame = 2)


# # Stop point #1 : please check the 'physical' coordinates in the images!
# import sys
# sys.exit('Image check ...')


# Mask coordinates
import copy


# ----- Format ----- #
#fix_sci{n} = {'slit1' : [[centers], [widths]],
#              'slit2' : [[centers], [widths]]}

# # 700 nm
# fix_sci01 = {'slit1' : [[2355, 1483, 1783, 1815],
#                         [28, 14, 26, 26, 24, 26]],
#              'slit2' : [[1750, 2310, 1800, 1837],
#                         [30, 13, 26, 26]]}
# fix_sci02 = {'slit1' : [[1482, 2355, 1783, 1815],
#                         [12, 28, 26, 26, 24, 26]],
#              'slit2' : [[1641, 1747, 2310, 1800, 1837],
#                         [50, 40, 14, 26, 26]]}

# 680 nm
fix_sci01 = {'slit1' : [[1483, 2355],
                        [12, 28]],
             'slit2' : [[1632, 1740, 2305],
                        [40, 50, 14]]}
fix_sci02 = {'slit1' : [[1484, 2355],
                        [14, 30]],
             'slit2' : [[492, 2305, 1633, 1741],
                        [12, 14, 50, 40]]}

# ------------------ #


for i in np.arange(len(iraf.type(ic.lst_sci, Stdout=1))):
    sci = iraf.type(ic.lst_sci, Stdout=1)[i].strip()
    exec('dic = copy.deepcopy(fix_sci{0:02d})'.format(i+1))
    for j in np.arange(ic.nslit):
    	fmsk = 'mskbadcol_'+sci+'_{0:1d}.txt'.format(j+1)
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

iraf.imdelete('xeqxbrg@'+ic.lst_sci)
iraf.imdelete('tmp@'+ic.lst_sci)
iraf.imdelete('tmpdq@'+ic.lst_sci)

for sci in iraf.type(ic.lst_sci, Stdout=1):
    sci = sci.strip()
    iraf.copy('eqxbrg'+sci+'.fits', 'tmp'+sci+'.fits')
    for j in np.arange(ic.nslit):
    	iraf.proto.fixpix('tmp'+sci+'.fits[sci,{0:1d}]'.format(j+1),
    		              'mskbadcol_'+sci+'_{0:1d}.pl'.format(j+1),
    		              linterp='1,2,3,4')
    iraf.copy('tmp'+sci+'.fits', 'xeqxbrg'+sci+'.fits')
    for j in np.arange(ic.nslit):
    	iraf.imarith('mskbadcol_'+sci+'_{0:1d}.pl'.format(j+1), '+',
    		         'xeqxbrg'+sci+'.fits[dq,{0:1d}]'.format(j+1),
    		         'tmpdq'+sci+'.fits[dq,{0:1d}]'.format(j+1))
    	iraf.imcopy('tmpdq'+sci+'.fits[dq,{0:1d}][*,*]'.format(j+1),
    		        'xeqxbrg'+sci+'.fits[dq,{0:1d}][*,*]'.format(j+1))
    	iraf.proto.fixpix('tmp'+sci+'.fits[sci,{0:1d}]'.format(j+1),
    		              'mskbadcol_'+sci+'_{0:1d}.pl'.format(j+1),
    		              linterp='1,2,3,4')    


for sci in iraf.type(ic.lst_sci, Stdout=1):
	for j in np.arange(ic.nslit):
		os.system('ds9 &')
		iraf.sleep(5.0)
		iraf.display(image = 'eqxbrg'+sci+'.fits[sci,{0:1d}]'.format(j+1), frame=1)
		iraf.display(image = 'xeqxbrg'+sci+'.fits[sci,{0:1d}]'.format(j+1), frame = 2)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))