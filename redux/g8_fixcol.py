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

        # 2-slit mode
        dt1, hd1 = fits.getdata('eqxbrg'+sci0+'.fits', ext=2, header=True)
        if (ic.nslit == 2):
            dt2, hd2 = fits.getdata('eqxbrg'+sci0+'.fits', ext=5, header=True)

        # Reading masking coordinates
        if (ic.nslit == 1):
            fix_sci = {'slit1': [[], []]}
        if (ic.nslit == 2):
            fix_sci = {'slit1': [[], []], 'slit2': [[], []]}

        for i in 1+np.arange(ic.nslit):
            regfile = "badcol/sci1_slit{0:1d}.reg".format(i)
            if (glob.glob(regfile) == []):
                fix_sci["slit{0:1d}".format(i)][0].append(100)
                fix_sci["slit{0:1d}".format(i)][1].append(0)
            else:
                x0, y0, xs, ys, _ = read_region(regfile, regtype="box")
                for j in np.arange(len(x0)):
                    fix_sci["slit{0:1d}".format(i)][0].append(int(round(x0[j])))
                    fix_sci["slit{0:1d}".format(i)][1].append(int(round(xs[j])))

        # Saving masking coordinates in text files
        dic = copy.deepcopy(fix_sci)
        for j in np.arange(ic.nslit):
            fmsk = 'mskbadcol_'+sci0+'_{0:1d}.txt'.format(j+1)
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

        # Masking bad columns with IRAF/fixpix
        iraf.copy('eqxbrg'+sci0+'.fits', 'tmp'+sci0+'.fits')
        for j in np.arange(ic.nslit):
            iraf.proto.fixpix('tmp'+sci0+'.fits[sci,{0:1d}]'.format(j+1),
                              'mskbadcol_'+sci0+'_{0:1d}.pl'.format(j+1),
                              linterp='1,2,3,4')
        iraf.copy('tmp'+sci0+'.fits', 'xeqxbrg'+sci0+'.fits')
        for j in np.arange(ic.nslit):
            iraf.imarith('mskbadcol_'+sci0+'_{0:1d}.pl'.format(j+1), '+',
                         'xeqxbrg'+sci0+'.fits[dq,{0:1d}]'.format(j+1),
                         'tmpdq'+sci0+'.fits[dq,{0:1d}]'.format(j+1))
            iraf.imcopy('tmpdq'+sci0+'.fits[dq,{0:1d}][*,*]'.format(j+1),
                        'xeqxbrg'+sci0+'.fits[dq,{0:1d}][*,*]'.format(j+1))
            iraf.proto.fixpix('tmp'+sci0+'.fits[sci,{0:1d}]'.format(j+1),
                              'mskbadcol_'+sci0+'_{0:1d}.pl'.format(j+1),
                              linterp='1,2,3,4') 

        # Checking if the masking task is well done
        z1l, z1u = np.percentile(dt1, [15, 85])
        # 2-slit mode
        if (ic.nslit == 2):
            z2l, z2u = np.percentile(dt2, [15, 85])

        if (ic.nslit == 1):
            z1, z2 = z1l, z1u
            ds9_frm = "ds9 xeqxbrg"+sci0+".fits[2] -multiframe"
            ds9_loc = " -scale lock yes -frame lock image"
            ds9_scl = " -scale limits {0:.2f} {1:.2f} &".format(z1, z2)
        if (ic.nslit == 2):
            z1, z2 = 0.5*(z1l+z2l), 0.5*(z1u+z2u)
            ds9_frm = "ds9 xeqxbrg"+sci0+".fits[2] xeqxbrg"+sci0+".fits[5] -multiframe"
            ds9_loc = " -scale lock yes -frame lock image"
            ds9_scl = " -scale limits {0:.2f} {1:.2f} &".format(z1, z2)

        os.system(ds9_frm + ds9_loc + ds9_scl)

        # Coming back to current path
        os.chdir(current_dir)
        iraf.chdir(current_dir)         


# Printing the running time
print('--- %.4f seconds ---' %(time.time()-start_time))
