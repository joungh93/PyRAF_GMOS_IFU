#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 17:21:30 2019

@author: jlee
"""


import numpy as np
import glob, os


# ----- Directory & file name & basic setting ----- #
procbias = 'Mbias.fits'
rawdir = '../raw/'
caldir = '../calibrations/'
dir_iraf = '../'
dir_db = 'database/'

lst_std = 'std.lis'
lst_arc = 'std_arc.lis'
lst_flat = 'std_flat.lis'

# Slit mode of IFU
nslit = 2


# ----- Find star w/ iraf.dir('onedstds') ----- #
## Absolute path: ~/[anaconda home directory]/envs/iraf27/iraf/noao/lib/onedstds/
## Check the standard star name from the header of the raw image
## $ find . -name *[standard starname]*


starname = 'gd108'
'''
Standard star name
exact name of [starname].dat
'''

stardir = 'onedstds$oke1990/'
'''
Directory name of standard star data file
'onedstds$[subdirectory]'
'''

extinction = 'onedstds$ctioextinct.dat'
'''
Extinction file name
GMOS-N: 'mk_extinct.txt' (needed to be downloaded)
GMOS-S: 'onedstds$ctioextinct.dat'
'''

root_name = starname+'_700_20190304_'
'''
Output file name for sensitivity function
(i.e. [starname]_[centwave]_[obsdate]_)
'''

obs_site = 'Gemini-South'
'''
Observing site
'Gemini-North' or 'Gemini-South'
'''


# ----- Check the mdf name w/ iraf.dir('gmos$data/*ifu*.fits', ncols=1) ----- #
## Absolute path: ~/[anaconda home directory]/envs/geminiconda/iraf_extern/gemini/gmos/data/
mdf = 'gsifu_slits_mdf_HAM.fits'
'''
Default MDF name
g[n/s]ifu_[ns]_slit[b/r/s]_mdf_[CCD].fits
n/s: Gemini north or south
[ns]: Nod & shuffle mode or not
[b/r/s]: blue/red/two slit mode
[CCD]: EEV or HAMAMATSU
'''

nmdf = 'new_gsifu_slits_mdf_HAM.fits'
'''
New MDF name
'''
