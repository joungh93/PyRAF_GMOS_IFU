# ----- Imports ----- #
import numpy as np
import glob, os
import pandas as pd
from pathlib import Path

from datetime import datetime
format_date_string = "%Y-%m-%d"

import warnings



#################################
##### User-defined variable #####
#################################

# Preferred number of bias frames
BIAS_TARGET_N = 30

# Search progressively farther from the science observation
BIAS_MARGINS = [3, 4, 5, 7, 10, 20]

#################################
#################################



# ----- Directories ----- #
dir_bias = Path('bias')
dir_std  = Path('standard')
dir_red  = Path('redux')

for di in [dir_bias, dir_std, dir_red]:
    di.mkdir(parents=True, exist_ok=True)
    
current_dir = Path.cwd().resolve()


# ----- Reading 'info.csv' ----- #
df = pd.read_csv("info.csv", sep="|")

dlab_id = df['DATALAB'].str[:-4].values
seq_num = df['DATALAB'].str[-3:].values.astype('int')


# ----- Data classification ----- #
bias = (df['OBJTYPE'].values == 'BIAS') & np.isin(df['OBSCLASS'].values, ["dayCal", "progCal"])
arc  = (df['OBJTYPE'].values == 'ARC')
flat = (df['OBJTYPE'].values == 'FLAT')
obj  = (df['OBJTYPE'].values == 'OBJECT')
sci  = (obj & (df['OBSCLASS'].values == 'science'))
std  = (obj & (df['OBSCLASS'].values == 'partnerCal'))


# ----- OBJECT & FLAT & ARC list ----- #
cen_wav, cnt_wav = np.unique(df['CENTWAVE'][sci].values, return_counts=True)
for i in np.arange(len(cen_wav)):
    dir_wav = Path(f"w{10*cen_wav[i]:.0f}")
    sci_wav = (sci & (df['CENTWAVE'].values == cen_wav[i]))
    sci_idx = np.where(sci_wav)[0]
    
    (dir_red / dir_wav).mkdir(parents=True, exist_ok=True)
    
    for j in np.arange(cnt_wav[i]):
        dir_sci = Path(df['FILENAME'][sci_wav].values[j])
        (dir_red / dir_wav / dir_sci).mkdir(parents=True, exist_ok=True)

        f = open(dir_red / dir_wav / dir_sci / 'sci.lis', 'w')
        f.write(f"{df['FILENAME'][sci_wav].values[j]}\n")
        f.close()

        f = open(dir_red / dir_wav / dir_sci / 'flat.lis', 'w')
        flat_wav = (flat & (df['CENTWAVE'].values == cen_wav[i]))
        flat_idx = np.abs(df['MJD-OBS'][flat_wav].values-df['MJD-OBS'][sci_wav].values[j]).argmin()
        f.write(f"{df['FILENAME'][flat_wav].values[flat_idx]}\n")
        f.close()        

        f = open(dir_red / dir_wav / dir_sci / 'arc.lis', 'w')
        arc_wav = (arc & (df['CENTWAVE'].values == cen_wav[i]))
        arc_idx = np.abs(df['MJD-OBS'][arc_wav].values-df['MJD-OBS'][sci_wav].values[j]).argmin()
        f.write(f"{df['FILENAME'][arc_wav].values[arc_idx]}\n")
        f.close()


# ----- STARNDARD & FLAT & ARC list ----- #
n_std = np.sum(std)
std_idx = np.where(std)[0]

if n_std == 0:
    raise ValueError("No standard star frame found!")

elif n_std > 1:
    print(f"WARNING: {n_std} standard star frames found.")
    print("All standard frames will be processed independently.")

std_cen_wav, std_cnt_wav = np.unique(df['CENTWAVE'][std].values, return_counts=True)
for i in np.arange(len(std_cen_wav)):
    dir_wav = Path(f"w{10*std_cen_wav[i]:.0f}")
    std_wav = (std & (df['CENTWAVE'].values == std_cen_wav[i]))
    std_idx = np.where(std_wav)[0]

    (dir_std / dir_wav).mkdir(parents=True, exist_ok=True)
    
    for j in np.arange(std_cnt_wav[i]):
        dir_frame = Path(df['FILENAME'][std_wav].values[j])
        (dir_std / dir_wav / dir_frame).mkdir(parents=True, exist_ok=True)

        f = open(dir_std / dir_wav / dir_frame / 'std.lis', 'w')
        f.write(f"{df['FILENAME'][std_wav].values[j]}\n")
        f.close()

        f = open(dir_std / dir_wav / dir_frame / 'std_flat.lis', 'w')
        std_flat_wav = (flat & (df['CENTWAVE'].values == std_cen_wav[i]))
        std_flat_idx = np.abs(df['MJD-OBS'][std_flat_wav].values-df['MJD-OBS'][std_idx].values[j]).argmin()
        f.write(f"{df['FILENAME'][std_flat_wav].values[std_flat_idx]}\n")
        f.close()        

        f = open(dir_std / dir_wav / dir_frame / 'std_arc.lis', 'w')
        std_arc_wav = (arc & (df['CENTWAVE'].values == std_cen_wav[i]))
        std_arc_idx = np.abs(df['MJD-OBS'][std_arc_wav].values-df['MJD-OBS'][std_idx].values[j]).argmin()
        f.write(f"{df['FILENAME'][std_arc_wav].values[std_arc_idx]}\n")
        f.close()


# ----- Grouping the BIAS frames ----- #

### Combine science & standard frames for assigning the bias epochs
bias_idx = np.where(bias)[0]
use_science_idx = np.where(sci | std)[0]
use_science_idx = use_science_idx[
    np.argsort(df['MJD-OBS'].values[use_science_idx])
]

### User functions
def group_into_epochs(indices, mjd, max_gap=3.0, max_span=7.0):
    """
    Dividing the epochs of observational data.
    """
    indices = np.asarray(indices)
    indices = indices[np.argsort(mjd[indices])]

    epochs = []

    current = [indices[0]]
    for idx in indices[1:]:
        previous_idx = current[-1]
        first_idx = current[0]

        gap  = (mjd[idx] - mjd[previous_idx])
        span = (mjd[idx] - mjd[first_idx])

        if (gap <= max_gap and span <= max_span):
            current.append(idx)

        else:
            epochs.append(np.asarray(current))
            current = [idx]

    epochs.append(np.asarray(current))

    return epochs


def select_bias_for_epoch(epoch_idx, bias_idx, mjd, date):
    """
    Selecting the bias frames for each epoch.
    """
    tmin = np.min(mjd[epoch_idx])
    tmax = np.max(mjd[epoch_idx])

    selected = None
    used_margin = None

    for margin in BIAS_MARGINS:
        good = ((mjd[bias_idx] >= tmin - margin) & \
                (mjd[bias_idx] <= tmax + margin))
        
        date_good = np.unique(date[bias_idx][good])
        date_cnd  = np.isin(date[bias_idx], date_good)

        candidate = bias_idx[good | date_cnd]

        if len(candidate) > 0:
            selected = candidate
            used_margin = margin

        if len(candidate) >= BIAS_TARGET_N:
            break

    return (selected, used_margin)


epochs = group_into_epochs(
    use_science_idx,
    df['MJD-OBS'].values,
    max_gap=3.0,
    max_span=7.0,
)

print("\n----- Epochs -----")
for i, idx in enumerate(epochs):
    print(i, df['DATE-OBS'].values[idx])

print("\n")

bias_groups = {}
for i, epoch_idx in enumerate(epochs):
    selected_bias, margin = \
        select_bias_for_epoch(
            epoch_idx,
            bias_idx,
            df['MJD-OBS'].values,
            df['DATE-OBS'].values,
        )

    print(
        f"Epoch {i+1:02d}:",
        f"{df['DATE-OBS'].values[epoch_idx][0]} - {df['DATE-OBS'].values[epoch_idx][-1]},",
        f"N_bias={len(selected_bias)},",
        f"margin={margin} days",
        f"(Bias: {df['DATE-OBS'].values[selected_bias][0]} - {df['DATE-OBS'].values[selected_bias][-1]})"
    )
    
    bias_groups[f"epoch{i+1:02d}"] = {}
    bias_groups[f"epoch{i+1:02d}"]['indices'] = selected_bias
    bias_groups[f"epoch{i+1:02d}"]['days_window']  = margin
    
    first_date = datetime.strptime(df['DATE-OBS'].values[selected_bias][0], format_date_string)
    last_date  = datetime.strptime(df['DATE-OBS'].values[selected_bias][-1], format_date_string)
    
    bias_groups[f"epoch{i+1:02d}"]['days_max_dt']  = (last_date - first_date).days

assert len(bias_groups.keys()) == len(epochs)


# ----- Bias mapping file for Science & Standard frames ----- #
with open(Path(dir_bias) / "master_bias_map.txt", "w") as mapfile:
    mapfile.write(
        "TARGET|DATE|MASTER_BIAS|LIST_FILE|"
        "NBIAS|WINDOW_DAY|MAX_DT_DAY\n"
    )

for i, epoch_idx in enumerate(epochs):
    epoch_tag = list(bias_groups.keys())[i]
    
    master_bias_name = f"Mbias_{epoch_tag}.fits"
    bias_list_name   = f"bias_{epoch_tag}.lis"
    
    with open(Path(dir_bias) / "master_bias_map.txt", "a") as mapfile:
        for _idx in epoch_idx:
            mapfile.write(
                f"{df['FILENAME'].values[_idx]}|"
                f"{df['DATE-OBS'].values[_idx]}|"
                f"{master_bias_name}|"
                f"{bias_list_name}|"
                f"{len(bias_groups[epoch_tag]['indices'])}|"
                f"{bias_groups[epoch_tag]['days_window']:d}|"
                f"{bias_groups[epoch_tag]['days_max_dt']:d}\n"
            )
        

# ----- BIAS list ----- #
for i, epoch_idx in enumerate(epochs):
    epoch_tag = list(bias_groups.keys())[i]
    
    bias_list_name = f"bias_{epoch_tag}.lis"
    bias_list_path = Path(dir_bias) / bias_list_name
    
    with open(bias_list_path, "w") as f:
        for _idx in epoch_idx:
            f.write(
                f"{df['FILENAME'].values[_idx]}\n"
            )
