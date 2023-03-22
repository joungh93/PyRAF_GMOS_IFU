# redux
(updated on 2023. 3. 22.)

## Description
Gemini GMOS/IFU reduction for science object data
* Science object pre-processing is similar to the workflow of the standard star processing.

## Workflow

### 0) Initial configurations
* Now move to each directory with science object data (`redux/w[CENT_WAVE]/[OBS_ID]/`), and check if all the list files (`sci.lis`, `flat.lis`, and `arc.lis`) have proper files. It is normal for each list file to include only one FITS file. Then, `g0_init_cfg.py` should be revised for the names of necessary files and the slit mode. _**For the detailed explanation for the revision, please refer to those comments in the code.**_

```
$ cd redux/
$ vi g0_init_cfg.py
(Revising)

$ ipython
(Python 2.7)
```

### 1) Verify the MDF
* Two Python codes are involved in this step. We have to revise `rev_mdf.py` manually to correct missing fibers. The comments that helps this process are written in the codes because there are a few :unamused: interactive tasks here.

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

* When running ``g1_vmdf.py``, we should check if the fibers are well assigned in the flat data. If there are some missing fibers that makes wrong MDF,  _**please take a note of the locations of the missing fibers**_. We have to revise them by revising & running the next code.

```
> run rev_mdf.py

Extracting slit 1
Recenter apertures for erg[FLAT]_1? ('yes')
Edit apertures for erg[FLAT]_1? ('yes')
(IRAF graphics displaying... please check the fibers visually.)
- "w" + "e" (left bottom) + "e" (right top) : zoom-in
- "w" + "a" : zoom-out
- "q" : quitting the interactive task

1) If unsatisfied,
Ctrl+C and revise this code again.

2) If satisfied with the MDF trace results,

Trace apertures for erg[FLAT]_1? ('yes')
Fit traced positions for erg[FLAT]_1 interactively? ('NO')
Write apertures for erg[FLAT]_1 to database ('yes')
Extract aperture spectra for erg[FLAT]_1? ('yes' (IFU-2 slit) / 'NO' (IFU-1 slit))
--> For the IFU-1 slit, this is the end.
Review extracted spectra from erg[FLAT]_1? ('NO')

Extracting slit 2
Find apertures for erg[FLAT]_2? ('yes')
(...repeating the tasks for slit 2...)
Extract aperture spectra for erg[FLAT]_1? ('NO')
--> 'GFEXTRACT exit status: error' message will appear, but it is not a fault.
```

### 2) Trace reference (without QE correction)
* Now the trace reference has to be extracted to define the position of the light on the detector. In the code, we have to check ``pk_line`` is identical throughout the workflow.

```
> run g2_trace_ref.py
```

### 3) Wavelength solution
* This step reduces arc data to produce the wavelength solution. A few interactive tasks are needed here, too.

```
> run g3_mk_wvsol.py

Examine identifications interactively? (Enter)
(IRAF graphics of spectrum displaying...)
"The spectrum window"
	- "w" + "e" (left bottom) + "e" (right top) : zoom-in
	- "w" + "a" : zoom-out
	- "d" :  delete the line
	- "m" : mark the line
	- "f" : jump to the parabola window
	- "q" : quitting the interactive task
"The parabola window"
	- "d" : jump to the spectrum window
	- "f" : fit the line again
	- "q" : return to the spectrum window
For the two-slit mode, you have to do the manual check twice.
Fit dispersion function interactively? (no|yes|NO|YES) ('NO'): Enter
Output file : erg[ARC].fits, database/aperg[ARC]_[1,2], database/iderg[ARC]_[1,2]
```

### 4) Reducing the lamp flat
* First, we have to find the gaps between the fiber bundles. Since the IRAF task ``gffindblocks`` is usually not working properly, we find them manually.

```
$ jupyter-notebook & (or jupyter-lab &)
(running the interactive tasks...)
```

* Next, we run the following code to reduce the flat data.
```
> run g4_proc_flat.py
```

### 5) Pre-processing of the science frames
* Using the master bias and the MDF we defined in the previous steps, we then pre-process the science frames of each science data. After running ``gfreduce``, we subtract the scattered light in the fiber gaps with ``gfscatsub``.
```
> run g5_preproc.py
```

### 6) Cosmic ray rejection
* We remove the cosmic ray of science frames using ``gemcrspec``.
```
> run g6_crrej.py
```

### 7) 

```
> run rv_w0.py

(IRAF graphics of spectrum displaying...)
"The spectrum window"
	- "l" : identify the night sky lines
	- "f" : fit for the mean velocity
	- "q" : write the parameters to the log file and quit the task
The header information (`CRVAL*`) would be revised after running this task.
```

:smiley_cat:❓ For more detailed instructions, please refer to the comments in the codes. :turtle::whale: 
