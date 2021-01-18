# PyRAF_GMOS_IFU
(updated on 2021. 1. 18.)

## Description
Gemini GMOS/IFU reduction & analysis package imported by PyRAF
* This package is only applicable to an observing program with a single field (not yet to multiple programs/fields... :crying_cat_face: :sweat_drops:).
* 

## Prerequisites
* Gemini GMOS/IFU raw data and the associated calibration data (bias) will be needed.
* The current versions of Python modules in this package are below.
  * ``numpy == 1.18.5``
  * ``pandas == 1.1.3``
  * ``astropy == 4.0.2``
* The following files should be in the working directory.
  * `login.cl` : PyRAF startup file
  * `lacos_spec.cl` : L.A.Cosmic task definition
* Before beginning, you have to check all the observation log files from [Gemini Data Archive](https://archive.gemini.edu/searchform). If there are some problematic data files (due to **_low counts_** or **_saturation_**), you should remove them from `./raw/` directory or move them to other paths.

## Subdirectories
* `./analysis/` : An analysis directory of the processed GMOS/IFU data
* `./bias/` : A reduction directory for bias data
* `./calibrations/` : A calibration data backup directory (for safety)
* `./obslog/` : Text files of observational log text files retrieved from [Gemini Data Archive](https://archive.gemini.edu/searchform)
* `./raw/` : Raw data from [Gemini Data Archive](https://archive.gemini.edu/searchform) and `obslog.py` from [GMOS Data Reduction Cookbook](http://ast.noao.edu/sites/default/files/GMOS_Cookbook/) with a few bugs revised
* `./redux/` : A reduction directory for object data
* `./standard/` : A reduction directory for standard star data

## Workflow
### 1) Initial data check using SQL
```
cd raw/
python obslog.py obsLog.sqlite3
sqlite3 obsLog.sqlite3

.table
.fullschema

# Check all the observation types
SELECT File,DateObs,Instrument,Object,ObsType,ObsClass,CcdBin,RoI,Disperser,CentWave,T_exp,use_me
FROM obslog GROUP BY ObsType;    # BIAS, ARC, FLAT, OBJECT

# Check all the bias frames
SELECT File,DateObs,Instrument,Object,ObsType,ObsClass,CcdBin,RoI,Disperser,CentWave,T_exp,use_me
FROM obslog WHERE ObsType='BIAS' GROUP BY File;

# Check all the arc frames
SELECT File,DateObs,Instrument,Object,ObsType,ObsClass,CcdBin,RoI,Disperser,CentWave,T_exp,use_me
FROM obslog WHERE ObsType='ARC' GROUP BY File;

# Check the GCAL flat frames (objects)
SELECT File,DateObs,Instrument,Object,ObsType,ObsClass,CcdBin,RoI,Disperser,CentWave,T_exp,use_me
FROM obslog WHERE ObsType='FLAT' GROUP BY File;

# Check all the object frames
SELECT File,DateObs,Instrument,Object,ObsType,ObsClass,CcdBin,RoI,Disperser,CentWave,T_exp,use_me
FROM obslog WHERE ObsType='OBJECT' AND ObsClass='science' GROUP BY File;

# Check the standard star frames
SELECT File,DateObs,Instrument,Object,ObsType,ObsClass,CcdBin,RoI,Disperser,CentWave,T_exp,use_me
FROM obslog WHERE ObsType='OBJECT' AND ObsClass='partnerCal' GROUP BY File;

# Check the GCAL flat frames (standards)
SELECT File,DateObs,Instrument,Object,ObsType,ObsClass,CcdBin,RoI,Disperser,CentWave,T_exp,use_me
FROM obslog WHERE ObsType='FLAT' GROUP BY File;

# Check the twilight flat frames
SELECT File,DateObs,Instrument,Object,ObsType,ObsClass,CcdBin,RoI,Disperser,CentWave,T_exp,use_me
FROM obslog WHERE ObsType='OBJECT' AND ObsClass='dayCal' GROUP BY File;
```

