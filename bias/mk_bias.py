# ----- Imports ----- #
import time
start_time = time.time()

import numpy as np
import glob, os
import shutil
import subprocess
from pathlib import Path


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
bias_lists = sorted(glob.glob("bias_*.lis"))
if len(bias_lists) == 0:
    raise RuntimeError("No bias_*.lis files found.")


# ----- Making master biases ----- #
for lst_bias in bias_lists:
    epoch_tag = Path(bias_lists[0]).stem.split("bias_")[1]
    procbias = (f"Mbias_{tag}.fits")


iraf.imdelete(procbias)
iraf.imdelete('g@'+lst_bias)
iraf.gbias('@'+lst_bias, procbias, rawpath=rawdir, fl_vardq='yes')
# iraf.copy(procbias, caldir)
os.system(f"cp -rpv {procbias} {caldir}")

os.system("rm -rfv tmp*")
iraf.imdelete('g@'+lst_bias)


# ----- Inspecting the processed bias ----- #
ds9_path = shutil.which('ds9')
subprocess.Popen([ds9_path])

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
