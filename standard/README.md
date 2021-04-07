# standard
(updated on 2021. 4. 7.)

## Description
Gemini GMOS/IFU reduction for standard star data
* Standard star reduction is needed for flux calibration of science observations.

## Workflow

### Initial configurations
Now move to the directory with standard star data, and check all the list files (`std*.lis`) have proper files. It is normal for each list file to include only one FITS file. Then, `g0_init_cfg.py` should be revised for the names of necessary files and the slit mode.

```
$ cd standard/
$ cat std[*].lis
$ vi g0_init_cfg.py
(Revising)

$ ipython
(Python 2.7)
```

### Verify the MDF
Two Python codes are involved in this step. You have to revise `rev_mdf.py` manually to correct missing fibers. The comments that helps this process are written in the codes because there are a few :unamused: interactive tasks here. _*For the detailed process, please refer to those comments in the codes.*_

```
> run g1_vmdf.py
> run rev_mdf.py
```

###


:smiley_cat:❓ For more detailed instructions, please refer to the comments in the codes. :turtle::whale: 
