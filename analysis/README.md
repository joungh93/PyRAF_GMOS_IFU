# Analysis
(updated on 2021. 9. 1.)

## Description
* Combining the reduced Gemini GMOS/IFU cubes
* Emission line analysis

## Workflow
* All the analysis processes are performed under Python 3 environment.
```
$ ipython
(Python 3.7)
```

### 1) Combining cubes
```
> run view_cube.py
> run check_offset.py
> run align_spat.py
> run align_spec.py
> run combine_cubes.py
```

### 2) Emission line analysis

  #### 2-1. Sky spectra
```
> run spec0_sky.py
```

#### 2-2. Blank-region subtraction
```
$ jupyter-notebook blank_init.ipynb &
```

#### 2-3. Initial plotting
```
> run plt0_collapsed.py
```
* To define the boundary between galaxy disks and tails, IRAF/ellipse task should be used here under Python 2 environment.
```
$ jupyter-notebook def_boundary.ipynb &
```

#### 2-4. Voronoi binning
```
$ jupyter-notebook vbin_init.ipynb &
> run spec0_vbin.py
```

#### 2-5. Emission line analysis
```
$ jupyter-notebook linefit_init.ipynb &
(Importing linefit.py)
> run spec1_linefit.py
$ jupyter-notebook contour_init.ipynb &
> run plt1_contour.py
> run plt1_linefit.py
```
