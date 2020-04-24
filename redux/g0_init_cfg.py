#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 17:21:30 2019

@author: jlee
"""


import numpy as np
import glob, os


# ----- File name & directory ----- #
procbias = 'Mbias.fits'
rawdir = '../Raw/'
caldir = '../calibrations/'
dir_iraf = '../'
dir_db = 'database/'

lst_sci = 'sci.lis'
lst_arc = 'arc.lis'
lst_flat = 'flat.lis'

nslit = 2

# For standard calibration
sensfunc = 'gd108_700_20190304_sens'
extinction = 'onedstds$ctioextinct.dat'
obs_site = 'Gemini-South'

# Check the mdf name w/ iraf.dir('gmos$data/*ifu*.fits', ncols=1)
mdf = 'gsifu_slits_mdf_HAM.fits'  
nmdf = 'new_gsifu_slits_mdf_HAM.fits'  
bias = '../calibrations/Mbias.fits'
