# ----- Imports ----- #
import numpy as np
import glob, os
from pathlib import Path
from astropy.io import fits
import pandas as pd
import subprocess, json


# ==============================
# Directories
# ==============================

current_dir = Path.cwd().resolve()
rootdir = current_dir.parent

dir_iraf = rootdir
rawdir   = rootdir / "raw"
caldir   = rootdir / "calibrations"
dir_db   = current_dir / "database"
dir_bias = rootdir / "bias"


# ==============================
# Input lists
# ==============================

dir_wav = []    # Central wavelength directory
for d in sorted(glob.glob("w*")):
	if os.path.isdir(d):
		dir_wav.append(d)
	else:
		raise ValueError("Please re-check the standard directories.")

lst_std  = "std.lis"
lst_arc  = "std_arc.lis"
lst_flat = "std_flat.lis"


# ==============================
# IFU configuration
# ==============================

# Reading the information
dict_frame = {}
for d in dir_wav:
    os.chdir(d)
    dict_frame[d] = []
    for dir_frame in sorted(glob.glob("*")):
        dict_frame[d].append(dir_frame)
    os.chdir(current_dir)

dfb = pd.read_csv(dir_bias / "master_bias_map.txt", sep="|")

# Slit mode of IFU
'''
Slit mode for IFU (as an input parameter of gfreduce)
IFU-1 slit: cslit = 'red' / 'blue'
IFU-2 slit: cslit = 'both'

Slit mode for IFU (as an input parameter of gfextract)
IFU-1 slit: eslit = 'red' / 'blue'
IFU-2 slit: eslit = '*'
'''

rawfile_ref = rawdir / f"{dict_frame[dir_wav[0]][0]}.fits"
header_ref  = fits.getheader(rawfile_ref, ext=0)

ifu_mask = header_ref.get("MASKNAME", "")
assert len(ifu_mask) > 0

ifu_mask_suffix = ifu_mask.split("IFU-")[1]

if (ifu_mask_suffix == "2"):
    nslit = 2
    cslit = "both"    # as an input parameter of gfreduce
    eslit = "*"    # as an input parameter of gfextract
elif (ifu_mask_suffix == "R"):
    nslit = 1
    cslit = "red"    # as an input parameter of gfreduce
    eslit = "red"    # as an input parameter of gfextract
elif (ifu_mask_suffix == "B"):
    nslit = 1
    cslit = "blue"    # as an input parameter of gfreduce
    eslit = "blue"    # as an input parameter of gfextract
else:
    raise ValueError("nslit must be 1 or 2.")


# # ==============================
# # Mask Definition File (MDF) configuration
# # ==============================

# pk_line = 1400    # Dispersion pixel used to identify fiber positions

# '''
# Check the mdf name w/ iraf.dir('gmos$data/*ifu*.fits', ncols=1) ----- #
# Absolute path: ~/[anaconda home directory]/envs/iraf27/iraf_extern/gemini/gmos/data/
# '''

# mdf = 'gnifu_slits_mdf.fits'
# '''
# Default MDF name
# g[n/s]ifu_[ns]_slit[b/r/s]_mdf_[CCD].fits
# n/s: Gemini north or south
# [ns]: Nod & shuffle mode or not
# [b/r/s]: blue/red/two slit mode
# [CCD]: EEV or HAMAMATSU
# '''

# nmdf = 'new_'+mdf
# '''
# New MDF name
# '''


# ============================================================
# Finding standard star in the IRAF database
# ============================================================
conda_results = subprocess.run(
    ['conda', 'env', 'list', '--json'],
    capture_output=True,
    text=True,
    check=True,
)
conda_data = json.loads(conda_results.stdout)
env_paths  = conda_data['envs']
env_path   = [Path(p) for p in env_paths if p.split("/")[-1] == "geminiconda"][0]


####################
####################
####################
####################
####################
####################


'''
Find star w/ iraf.dir('onedstds') ----- #
Absolute path: ~/[anaconda home directory]/envs/iraf27/iraf/noao/lib/onedstds/
Check the standard star name from the header of the raw image (extension: 0, keyword: 'OBJECT')
$ find . -name *[standard starname]*
'''

starname = header_ref['OBJECT'].lower()    # 'wolf1346'
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


