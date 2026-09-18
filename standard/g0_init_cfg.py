# ----- Imports ----- #
import numpy as np
import glob, os
from pathlib import Path
from astropy.io import fits
import pandas as pd
import subprocess, json


# ==============================
# User functions
# ==============================

def get_mdf_name(
    rawfile,
    gmos_data_dir="/home/jhlee/miniconda3/envs/geminiconda/iraf_extern/gemini/gmos/data",
):
    with fits.open(rawfile) as hdul:
        hdr = hdul[0].header

    # Gemini North / South
    instrume = hdr["INSTRUME"].upper()
    if instrume == "GMOS-N":
        site = "n"
    elif instrume == "GMOS-S":
        site = "s"
    else:
        raise ValueError(f"Unknown instrument: {instrume}")

    # Nod & Shuffle
    nodmode = str(hdr.get("NODMODE", "")).strip().upper()
    if nodmode == "STANDARD":
        ns = "ns_"
    else:
        ns = ""

    # IFU slit mode
    maskname = hdr["MASKNAME"].strip().upper()
    slit_map = {
        "IFU-B": "slitb",
        "IFU-B-NS": "slitb",
        "IFU-R": "slitr",
        "IFU-R-NS": "slitr",
        "IFU-2": "slits",
        "IFU-NS-2": "slits",
    }

    if maskname not in slit_map:
        raise ValueError(f"Unknown IFU MASKNAME: {maskname}")

    slit = slit_map[maskname]

    # Detector type
    dettype = str(
        hdr.get("DETECTOR", hdr.get("DETTYPE", ""))
    ).upper()

    # date_obs = str(
        # hdr.get("DATE-OBS", "")
    # ).strip()

    if "HAMAMATSU_NEW" in dettype or "HAMAMATSU NEW" in dettype:
        ccd = "HAM-2"

    elif "HAMAMATSU" in dettype:
        ccd = "HAM"

    elif (
        "EEV" in dettype
        or "E2V" in dettype
        or "SDSU" in dettype
    ):
        ccd = None

    else:
        raise ValueError(
            f"Unknown detector: {dettype}"
        )

    if ccd is None:
        mdf = f"g{site}ifu_{ns}{slit}_mdf.fits"
    else:
        mdf = f"g{site}ifu_{ns}{slit}_mdf_{ccd}.fits"

    if not os.path.exists(os.path.join(gmos_data_dir, mdf)):
        raise FileNotFoundError(f"MDF file is not found: {mdf}")
        
    return mdf
    
    
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

obs_site = header_ref.get("TELESCOP", "")
'''
Observing site
'Gemini-North' or 'Gemini-South'
'''
assert len(obs_site) > 0


# ==============================
# Mask Definition File (MDF) configuration
# ==============================

conda_results = subprocess.run(
    ['conda', 'env', 'list', '--json'],
    capture_output=True,
    text=True,
    check=True,
)
conda_data = json.loads(conda_results.stdout)
env_paths  = conda_data['envs']
env_path   = [Path(p) for p in env_paths if p.split("/")[-1] == "geminiconda"][0]

pk_line = 1400    # Dispersion pixel used to identify fiber positions

'''
Check the mdf name w/ iraf.dir('gmos$data/*ifu*.fits', ncols=1) ----- #
Absolute path: ~/[anaconda home directory]/envs/iraf27/iraf_extern/gemini/gmos/data/
'''

mdf = get_mdf_name(
    rawfile_ref,
    gmos_data_dir=str(env_path / "iraf_extern" / "gemini" / "gmos" / "data"))
'''
Default MDF name
g[n/s]ifu_[ns]_slit[b/r/s]_mdf_[CCD].fits
n/s: Gemini north or south
[ns]: Nod & shuffle mode or not
[b/r/s]: blue/red/two slit mode
[CCD]: EEV or HAMAMATSU
'''

# Working MDF used by all later reduction steps.
nmdf = "new_" + mdf

# MDF review parameters
mdf_bundle_size = 50

# Only used to highlight unusually large separations
# in the diagnostic plot.
# This does NOT automatically reject any fiber.
mdf_gap_factor = 1.6

# User-selected additional missing fibers
mdf_missing_file = "mdf_missing_fibers.txt"

# Review history
mdf_review_log = "mdf_review.log"

# Saved diagnostic figures
mdf_review_prefix = "mdf_review"

# Set True only when you deliberately want to start
# the MDF inspection again from the original Gemini MDF.
mdf_reset = False


# ==============================
# Finding standard star in the IRAF database
# ==============================

'''
Find star w/ iraf.dir('onedstds') ----- #
Absolute path: ~/[anaconda home directory]/envs/iraf27/iraf/noao/lib/onedstds/
Check the standard star name from the header of the raw image (extension: 0, keyword: 'OBJECT')
$ find . -name *[standard starname]*
'''

try:
    starname = header_ref['OBJECT'].lower()    # e.g., 'wolf1346'
    '''
    Standard star name
    exact name of [starname].dat
    '''
    find_results = subprocess.run(
        ["find", f"{env_path}/iraf/noao/lib/onedstds/", "-name", f"*{starname}*"],
        capture_output=True,
        text=True,
        check=True,
    )
    find_star_list = find_results.stdout.split("\n")
    find_star_list.remove("")

    if len(find_star_list) == 0:
        raise ValueError(
            "Please manually check the star directory:\n"
            f"Move to {env_path}/iraf/noao/lib/onedstds/"
        )
        
    elif len(find_star_list) == 1:
        idx_select = 0
        
    else:
        nline_star_list = []
        for file in find_star_list:
            with open(file, 'r') as f:
                ll = f.readlines()
                nline_star_list.append(len(ll))
        nline_star_list = np.asarray(nline_star_list)
        idx_select = np.argmax(nline_star_list)

    subdir  = find_star_list[idx_select].split("onedstds/")[1].split("/")[0]
    stardir = f"onedstds${subdir}"
    
except:
    starname = None
    stardir  = None    # 'onedstds$spec50cal/'
    print("'starname' and 'stardir' are manually set.")
    print("  Current setting:")
    print(f"    'starname': {starname}")
    print(f"    'stardir': {stardir}")
    '''
    Directory name of standard star data file
    'onedstds$[subdirectory]/'
    '''

# root_name = starname+'_700_20190613_'
# '''
# Output file name for sensitivity function
# (i.e. [starname]_[centwave]_[obsdate]_)
# '''

if (obs_site.split('-')[1] == "North"):
    extinction = dir_iraf+'mk_extinct.txt'
elif (obs_site.split('-')[1] == "South"):
    extinction = 'onedstds$ctioextinct.dat'
'''
Extinction file name
GMOS-N: 'mk_extinct.txt' (needed to be downloaded from Buton+13)
GMOS-S: 'onedstds$ctioextinct.dat'
'''


# ==============================
# Save configuration
# ==============================

# If needed, please revise the following configurations manually!
config = {
    "dir_iraf": str(dir_iraf),
    "rawdir": str(rawdir),
    "caldir": str(caldir),
    "dir_db": str(dir_db),
    "dir_bias": str(dir_bias),

    "dir_wav": dir_wav,
    "dict_frame": dict_frame,

    "lst_std": lst_std,
    "lst_arc": lst_arc,
    "lst_flat": lst_flat,

    "nslit": nslit,
    "cslit": cslit,
    "eslit": eslit,

    "obs_site": obs_site,

    "pk_line": pk_line,

    "mdf": mdf,
    "nmdf": nmdf,

    "mdf_bundle_size": mdf_bundle_size,
    "mdf_gap_factor": mdf_gap_factor,
    "mdf_missing_file": mdf_missing_file,
    "mdf_review_log": mdf_review_log,
    "mdf_review_prefix": mdf_review_prefix,
    "mdf_reset": mdf_reset,

    "starname": starname,
    "stardir": stardir,
    "extinction": str(extinction),
}

with open("config.json", "w") as f:
    json.dump(config, f, indent=4)

print("Configuration saved to config.json")
