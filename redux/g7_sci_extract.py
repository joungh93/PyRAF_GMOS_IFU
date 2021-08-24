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

iraf.unlearn('gqecorr')
iraf.unlearn('gfextract')
iraf.unlearn('gfdisplay')


# ---------- QE correction and extraction of the science ---------- #
for d in ic.dir_wav:
    dir_sci = sorted(glob.glob(d+"/*"))

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

        # Coming back to current path
        os.chdir(current_dir)
        iraf.chdir(current_dir) 



# flat0 = iraf.type(ic.lst_flat, Stdout=1)[0].strip()
# response = flat0+'_resp'
# ref_flat0 = 'eqbrg'+flat0

# arc0 = iraf.type(ic.lst_arc, Stdout=1)[0].strip()

# iraf.imdelete('qxbrg@'+ic.lst_sci)
# iraf.imdelete('eqxbrg@'+ic.lst_sci)

# for sci in iraf.type(ic.lst_sci, Stdout=1):
#     sci = sci.strip()
#     iraf.gqecorr('xbrg'+sci, refimage='erg'+arc0, fl_correct='yes',
#                  fl_vardq='yes', verbose='yes')
#     iraf.gfextract('qxbrg'+sci, response=response, recenter='no',
#                    trace='no', reference=ref_flat0, weights='none',
#                    fl_vardq='yes', line=1400)

# os.system('ds9 &')
# iraf.sleep(5.0)
# for sci in iraf.type(ic.lst_sci, Stdout=1):
# 	sci = sci.strip()
# 	iraf.gfdisplay('eqxbrg'+sci, 1)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
