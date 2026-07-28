#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 17:21:30 2019

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import glob, os
import subprocess


# ----- File name & directory ----- #
cpath = os.path.abspath(".")
dir_iraf = "/".join(cpath.split("/")[:-1])+"/"
rawdir = dir_iraf+"raw/"
caldir = dir_iraf+"calibrations/"
if (glob.glob(caldir) == []):
	os.mkdir(caldir)
lst_bias = 'bias.lis'
procbias = 'Mbias.fits'


# ----- Importing IRAF from the root directory ----- #
current_dir = os.getcwd()
os.chdir(dir_iraf)

from pyraf import iraf
from pyraf.iraf import gemini, gmos

iraf.unlearn()
os.chdir(current_dir)
iraf.chdir(current_dir)


# ----- Making a processed bias ----- #
iraf.imdelete(procbias)
iraf.imdelete('g@'+lst_bias)
iraf.gbias('@'+lst_bias, procbias, rawpath=rawdir, fl_vardq='yes')
# iraf.copy(procbias, caldir)
os.system(f"cp -rpv {procbias} {caldir}")


# ----- Inspecting the processed bias ----- #
# os.system('/home/jhlee/Downloads/ds9 &')
subprocess.Popen(['/home/jhlee/Downloads/ds9'])

# Wait until DS9 is ready
for i in range(30):
    try:
        result = subprocess.run(
            ['xpaget', 'ds9',],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=1
        )

        if result.returncode == 0:
            break

    except Exception:
        pass

    time.sleep(1)

else:
    raise RuntimeError("DS9 did not become ready within 30 seconds.")

# iraf.sleep(5.0)
iraf.gdisplay(procbias, 1, fl_paste='no')


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))
