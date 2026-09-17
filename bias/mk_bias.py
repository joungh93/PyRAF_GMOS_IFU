# ----- Imports ----- #
import time
start_time = time.time()

import numpy as np
import glob, os
import shutil
import subprocess
from pathlib import Path
import pandas as pd


# ----- File name & directory ----- #
current_dir = Path.cwd().resolve()
dir_iraf = Path("../").resolve()
rawdir = dir_iraf / "raw"
caldir = dir_iraf / "calibrations"
caldir.mkdir(parents=True, exist_ok=True)

lst_bias = 'bias.lis'
procbias = 'Mbias.fits'


# ----- Importing IRAF from the root directory ----- #
os.chdir(dir_iraf)

from pyraf import iraf
from pyraf.iraf import gemini, gmos

os.chdir(current_dir)
iraf.chdir(f"{current_dir}")

iraf.unlearn("gbias")


# ----- Finding all bias lists ----- #
dfb = pd.read_csv("master_bias_map.txt", sep="|")

bias_lists = sorted(glob.glob("bias_*.lis"))
if len(bias_lists) == 0:
    raise RuntimeError("No bias_*.lis files found.")

assert bias_lists == list(np.unique(dfb['LIST_FILE'].values))


# ----- Making master biases ----- #
for lst_bias in bias_lists:
    epoch_tag = Path(lst_bias).stem.split("bias_")[1]
    procbias  = f"Mbias_{epoch_tag}.fits"
    
    with open(lst_bias, 'r') as f:
        bias_names = [line.strip() for line in f if line.strip()]
    nbias = len(bias_names)
    
    print()
    print(f"--- Making {procbias} ---")
    print(f"--- Number of bias frames: {nbias} ---")
    
    if nbias < 30:
        print(
            "WARNING: fewer than 30 bias frames "
            "are being combined."
        )

    if nbias < 10:
        print(
            "WARNING: very small bias sample."
        )
    
    ### Remove old output if present
    if os.path.exists(procbias):
        iraf.imdelete(procbias, verify="no")
    
    ### Remove temporary/prepared bias files
    try:
        iraf.imdelete("g@"+lst_bias, verify="no")
    except Exception:
        pass
    
    ### Run gbias task
    logfile = f"gbias_{epoch_tag}.log"
    
    iraf.gbias(
        "@"+lst_bias,
        procbias,
        rawpath=f"{rawdir}",
        logfile=logfile,
    )
    
    ### Copy to calibration directory
    shutil.copy2(procbias, caldir / procbias)
    
    ### Clean only temporary products
    os.system("rm -rfv tmp*")
    try:
        iraf.imdelete("g@"+lst_bias, verify="no")
    except Exception:
        pass
        

# # ----- Inspecting the processed bias ----- #
# ds9_path = shutil.which('ds9')
# subprocess.Popen([ds9_path])

# # Wait until DS9 is ready
# for i in range(30):
    # try:
        # result = subprocess.run(
            # ['xpaget', 'ds9',],
            # stdout=subprocess.PIPE,
            # stderr=subprocess.PIPE,
            # timeout=1
        # )

        # if result.returncode == 0:
            # break

    # except Exception:
        # pass

    # time.sleep(1)

# else:
    # raise RuntimeError("DS9 did not become ready within 30 seconds.")

# # iraf.sleep(5.0)
# iraf.gdisplay(procbias, 1, fl_paste='no')


# ----- Print the running time ----- #
print()
print(
    "--- Total running time: "
    f"{time.time() - start_time:.1f} sec ---"
)
