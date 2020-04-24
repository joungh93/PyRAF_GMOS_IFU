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

iraf.unlearn('gemcrspec')


# ---------- Cosmic ray rejection ---------- #
iraf.imdelete('xbrg@'+ic.lst_sci)

for sci in iraf.type(ic.lst_sci, Stdout=1):
    sci = sci.strip()
    iraf.gemcrspec('brg'+sci, 'xbrg'+sci, logfile='crrej.log',
                   key_gain='GAIN', key_ron='RDNOISE', xorder=9,
                   yorder=-1, sigclip=4.5, sigfrac=0.5, objlim=1.,
                   niter=4, verbose='yes', fl_vardq='yes')

# os.system('ds9 &')
# iraf.sleep(5.0)
# for sci in iraf.type(ic.lst_sci, Stdout=1):
#     sci = sci.strip()
#     for i in range(12):
#         iraf.imexamine('xbrg'+sci+'[sci,'+str(i+1)+']', 1)


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
