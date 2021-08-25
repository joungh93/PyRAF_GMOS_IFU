#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 3 01:18:27 2020

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


# ----- Wavelength alignment ----- #
from astropy.io import fits
import pandas as pd

iraf.rv()
iraf.unlearn('rvidlines')
iraf.rv.observatory = ic.obs_site

for d in ic.dir_wav:
    dir_sci = sorted(glob.glob(d+"/*"))

    for j in np.arange(len(dir_sci)):

        # Moving each science directory
        name_sci = dir_sci[j].split("/")[-1]
        print("Moving path for "+name_sci+"...")
        os.chdir(current_dir+"/"+dir_sci[j])
        iraf.chdir(current_dir+"/"+dir_sci[j])

        # SCI
        sci = np.loadtxt(ic.lst_sci, dtype=str)
        sci0 = sci.item(0)

        # Running rvidlines
        hd0 = fits.getheader('stxeqxbrg'+sci0+'.fits', ext=0)    # Primary HDU
        hd2 = fits.getheader('stxeqxbrg'+sci0+'.fits', ext=2)    # SCI HDU

        os.system('rm -rfv rvLog_'+sci0+'.txt')              
        iraf.rvidlines('stxeqxbrg'+sci0+'.fits[SKY]', coordlist=ic.dir_iraf+'skylines.txt', ftype='emission',
                       nsum=1, threshold=7.0, maxfeatures=13, fwidth=10.0,
                       cradius=10.0, minsep=5.0, autowrite='yes',
                       logfile='rvLog_'+sci0+'.txt')

        f = open('rvLog_'+sci0+'.txt','r')
        lines = f.readlines()
        bool_Zobs = pd.Series(lines).str.startswith('stxeqxbrg'+sci0+'.fits[SKY]   1 : Zobs').values
        line_Zobs = np.array(lines)[bool_Zobs][0]
        Zobs = float(line_Zobs.split('Zobs     =')[1].split(',')[0])
        f.close()

        print('{0:.3f} --> {1:.3f}'.format(hd2['CRVAL1'], hd2['CRVAL1']+(-Zobs*hd0['CENTWAVE']*10.0)))
        nw_zero = hd2['CRVAL1']+(-Zobs*hd0['CENTWAVE']*10.0)

        iraf.unlearn('hedit')
        iraf.hedit.update = 'yes'
        iraf.hedit.verify = 'no'
        iraf.hedit('stxeqxbrg'+sci0+'.fits[SCI]', 'CRVAL1', nw_zero)

        # Coming back to current path
        os.chdir(current_dir)
        iraf.chdir(current_dir)  
