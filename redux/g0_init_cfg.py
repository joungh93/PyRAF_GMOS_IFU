#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 17:21:30 2019

@author: jlee
"""


import numpy as np
import glob, os


# ----- Directory & file name & basic setting ----- #
cpath = os.path.abspath(".")
dir_iraf = "/".join(cpath.split("/")[:-1])+"/"
rawdir = dir_iraf+"raw/"
caldir = dir_iraf+"calibrations/"
dir_db = "database/"    # Relative path
dir_std = dir_iraf+"standard/"
dir_wav = []    # Central wavelength directory
for d in sorted(glob.glob("w*")):
	if os.path.isdir(d):
		dir_wav.append(d)
	else:
		raise ValueError("Please re-check the science directories.")

procbias = 'Mbias.fits'

lst_std = 'std.lis'
lst_stdarc = 'std_arc.lis'
lst_stdflat = 'std_flat.lis'

lst_sci = 'sci.lis'
lst_arc = 'arc.lis'
lst_flat = 'flat.lis'
pk_line = 1400    # Deteremined from standard star reduction
'''
Line for finding peaks (gfreduce)
Line/column for finding apertures (gfextract)
'''

# Slit mode of IFU
nslit = 2
cslit = 'both'
'''
Slit mode for IFU (as an input parameter of gfreduce)
IFU-1 slit: cslit = 'red' / 'blue'
IFU-2 slit: cslit = 'both'
'''

if (nslit == 1):
	eslit = cslit
if (nslit == 2):
	eslit = '*'
'''
Slit mode for IFU (as an input parameter of gfextract)
IFU-1 slit: eslit = 'red' / 'blue'
IFU-2 slit: eslit = '*'
'''


# ----- Find star w/ iraf.dir('onedstds') ----- #
## Absolute path: ~/[anaconda home directory]/envs/iraf27/iraf/noao/lib/onedstds/
## Check the standard star name from the header of the raw image (extension: 0, keyword: 'OBJECT')
## $ find . -name *[standard starname]*


starname = 'feige34'
'''
Standard star name
exact name of [starname].dat
'''

stardir = 'onedstds$spec50cal/'
'''
Directory name of standard star data file
'onedstds$[subdirectory]/'
'''

extinction = dir_iraf+'mk_extinct.txt'
'''
Extinction file name
GMOS-N: 'mk_extinct.txt' (needed to be downloaded from Buton+13)
GMOS-S: 'onedstds$ctioextinct.dat'
'''

sensfunc = starname+'_650_20210408_sens'
'''
Output file name for sensitivity function
(i.e. [starname]_[centwave]_[obsdate]_sens)
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
