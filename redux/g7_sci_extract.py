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
from astropy.io import fits


# ----- Importing IRAF from the root directory ----- #
current_dir = os.getcwd()
os.chdir(ic.dir_iraf)

from pyraf import iraf
from pyraf.iraf import gemini, gmos

os.chdir(current_dir)
iraf.chdir(current_dir)

iraf.unlearn('gqecorr')
iraf.unlearn('gfextract')
iraf.unlearn('gfdisplay')


# ---------- QE correction and extraction of the science ---------- #
for d in ic.dir_wav:
    dir_sci = sorted(glob.glob(d+"/*"))

    filenames1 = ""
    if (ic.nslit == 2):
        filenames2 = ""

    for j in np.arange(len(dir_sci)):

        # Moving each science directory
        name_sci = dir_sci[j].split("/")[-1]
        print("Moving path for "+name_sci+"...")
        os.chdir(current_dir+"/"+dir_sci[j])
        iraf.chdir(current_dir+"/"+dir_sci[j])

        # SCI & ARC & FLAT
        sci = np.loadtxt(ic.lst_sci, dtype=str)
        sci0 = sci.item(0)

        arc = np.loadtxt(ic.lst_arc, dtype=str)
        arc0 = arc.item(0)

        flat = np.loadtxt(ic.lst_flat, dtype=str)
        flat0 = flat.item(0)
        response = flat0+'_resp'

        # QE correction
        iraf.imdelete('qxbrg@'+ic.lst_sci)
        iraf.imdelete('eqxbrg@'+ic.lst_sci)
        iraf.gqecorr('xbrg'+sci0, refimage='erg'+arc0, fl_correct='yes',
                     fl_vardq='yes', verbose='yes')

        # Science extraction
        iraf.gfextract('qxbrg'+sci0, response=response, recenter='no',
                       trace='no', reference='eqbrg'+flat0, weights='none',
                       fl_vardq='yes', line=ic.pk_line)        

        # Displaying the extracted data for checking bad columns
        if (glob.glob("badcol/") == []):
            os.system("mkdir badcol")
        # ---> making a directory for bad column masking region files
        # ("badcol/sci{}_slit{}.reg", with DS9 saoimage format)

        filenames1 += (dir_sci[j]+"/eqxbrg"+sci0+".fits[2] ")
        if (ic.nslit == 2):
            filenames2 += (dir_sci[j]+"/eqxbrg"+sci0+".fits[5] ")

        # #####
        # dt1, hd1 = fits.getdata('eqxbrg'+sci0+'.fits', ext=2, header=True)
        # z1l, z1u = np.percentile(dt1, [15, 85])
        # # 2-slit mode
        # if (ic.nslit == 2):
        #     dt2, hd2 = fits.getdata('eqxbrg'+sci0+'.fits', ext=5, header=True)
        #     z2l, z2u = np.percentile(dt2, [15, 85])

        # if (ic.nslit == 1):
        #     z1, z2 = z1l, z1u
        #     ds9_frm = "ds9 eqxbrg"+sci0+".fits[2] -multiframe"
        #     ds9_loc = " -scale lock yes -frame lock image"
        #     ds9_scl = " -scale limits {0:.2f} {1:.2f} &".format(z1, z2)
        # if (ic.nslit == 2):
        #     z1, z2 = 0.5*(z1l+z2l), 0.5*(z1u+z2u)
        #     ds9_frm = "ds9 eqxbrg"+sci0+".fits[2] eqxbrg"+sci0+".fits[5] -multiframe"
        #     ds9_loc = " -scale lock yes -frame lock image"
        #     ds9_scl = " -scale limits {0:.2f} {1:.2f} &".format(z1, z2)

        # os.system(ds9_frm + ds9_loc + ds9_scl)
        # #####

        # Coming back to current path
        os.chdir(current_dir)
        iraf.chdir(current_dir)

    ds9_opt = "-scalemode zscale -scale lock yes -frame lock image"
    os.system("ds9 "+ds9_opt+" "+filenames1+"&")
    if (ic.nslit == 2):
        os.system("ds9 "+ds9_opt+" "+filenames2+"&")


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
