# standard
(updated on 2021. 4. 12.)

## Description
Gemini GMOS/IFU reduction for standard star data
* Standard star reduction is needed for flux calibration of science observations.

## Workflow

### 1) Initial configurations
Now move to the directory with standard star data, and check all the list files (`std*.lis`) have proper files. It is normal for each list file to include only one FITS file. Then, `g0_init_cfg.py` should be revised for the names of necessary files and the slit mode. _**For the detailed explanation for the revision, please refer to those comments in the code.**_

```
$ cd standard/
$ cat std[*].lis
$ vi g0_init_cfg.py
(Revising)

$ ipython
(Python 2.7)
```

### 2) Verify the MDF
Two Python codes are involved in this step. You have to revise `rev_mdf.py` manually to correct missing fibers. The comments that helps this process are written in the codes because there are a few :unamused: interactive tasks here.

```
> run g1_vmdf.py

Extracting slit 1
Find apertures for erg[FLAT]_1? ('yes')
Edit apertures for erg[FLAT]_1? ('yes')
- "w" + "e" (left bottom) + "e" (right top) : zoom-in
- "w" + "a" : zoom-out
- "q" : quitting the interactive task
Trace apertures for erg[FLAT]_1? ('yes')
Fit traced positions for erg[FLAT]_1 interactively? ('NO')
Write apertures for erg[FLAT]_1 to database ('yes')
Extract aperture spectra for erg[FLAT]_1? ('yes' (IFU-2 slit) / 'NO' (IFU-1 slit))

Review extracted spectra from erg[FLAT]_1? ('NO')

Extracting slit 2
Find apertures for erg[FLAT]_2? ('yes')
Extract aperture spectra for erg[FLAT]_1? ('NO')
```

When running ``g1_vmdf.py``, you should check if the fibers are well assigned in the flat data. If there are some missing fibers that makes wrong MDF,  _**please take a note of the locations of the missing fibers**_. We have to revise them by revising & running the next code.

```
> run rev_mdf.py
```

###


:smiley_cat:❓ For more detailed instructions, please refer to the comments in the codes. :turtle::whale: 
