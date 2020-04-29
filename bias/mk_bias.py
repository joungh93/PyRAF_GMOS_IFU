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


# ----- File name & directory ----- #
lst_bias = 'bias.lis'
procbias = 'Mbias.fits'
rawdir = '../Raw/'
caldir = '../calibrations/'


# ----- Importing IRAF from the root directory ----- #
current_dir = os.getcwd()
os.chdir('../')

from pyraf import iraf
from pyraf.iraf import gemini, gmos

iraf.unlearn()
os.chdir(current_dir)
iraf.chdir(current_dir)


# ----- Making a processed bias ----- #
iraf.imdelete(procbias)
iraf.imdelete('g@'+lst_bias)
iraf.gbias('@'+lst_bias, procbias, rawpath=rawdir, fl_vardq='yes')
iraf.copy(procbias, caldir)


# ----- Inspecting the processed bias ----- #
os.system('ds9 &')
iraf.sleep(5.0)
iraf.gdisplay(procbias, 1, fl_paste='no')


# Printing the running time
print('--- %s seconds ---' %(time.time()-start_time))