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
cslit = 'both'
'''
Slit mode for IFU
IFU-2 slit: cslit = 'both'
IFU-1 slit: cslit = 'red' / 'blue'
'''


# ----- Find star w/ iraf.dir('onedstds') ----- #
## Absolute path: ~/[anaconda home directory]/envs/iraf27/iraf/noao/lib/onedstds/
## Check the standard star name from the header of the raw image (extension: 0, keyword: 'OBJECT')
## $ find . -name *[standard starname]*


starname = 'wolf1346'
'''
Standard star name
exact name of [starname].dat
'''

stardir = 'onedstds$spec50cal/'
'''
Directory name of standard star data file
'onedstds$[subdirectory]/'
'''

extinction = 'mk_extinct.txt'
'''
Extinction file name
GMOS-N: 'mk_extinct.txt' (needed to be downloaded)
GMOS-S: 'onedstds$ctioextinct.dat'
'''

root_name = starname+'_700_20190613_'
'''
Output file name for sensitivity function
(i.e. [starname]_[centwave]_[obsdate]_)
'''

obs_site = 'Gemini-North'
'''
Observing site
'Gemini-North' or 'Gemini-South'
'''


# ----- Check the mdf name w/ iraf.dir('gmos$data/*ifu*.fits', ncols=1) ----- #
## Absolute path: ~/[anaconda home directory]/envs/iraf27/iraf_extern/gemini/gmos/data/
mdf = 'gnifu_slits_mdf.fits'
'''
Default MDF name
g[n/s]ifu_[ns]_slit[b/r/s]_mdf_[CCD].fits
n/s: Gemini north or south
[ns]: Nod & shuffle mode or not
[b/r/s]: blue/red/two slit mode
[CCD]: EEV or HAMAMATSU
'''

nmdf = 'new_'+mdf
'''
New MDF name
'''
