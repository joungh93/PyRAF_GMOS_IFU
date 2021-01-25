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
rawdir = '../raw/'
caldir = '../calibrations/'
dir_iraf = '../'
dir_db = 'database/'

lst_std = 'std.lis'
lst_arc = 'std_arc.lis'
lst_flat = 'std_flat.lis'

nslit = 2


# Find star w/ iraf.dir('onedstds') or iraf.dir('onedstds$oke1990')
starname = 'gd108'
stardir = 'onedstds$oke1990/'
extinction = 'onedstds$ctioextinct.dat'
root_name = starname+'_700_20190304_'
obs_site = 'Gemini-South'

# Check the mdf name w/ iraf.dir('gmos$data/*ifu*.fits', ncols=1)
mdf = 'gsifu_slits_mdf_HAM.fits'
nmdf = 'new_gsifu_slits_mdf_HAM.fits'    
bias = '../calibrations/Mbias.fits'

