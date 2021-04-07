# bias
(updated on 2021. 4. 7.)

## Description
Gemini GMOS/IFU reduction for bias data

## Workflow

### Creating master bias
The only code for bias reduction is ``mk_bias.py``.

```
$ conda activiate iraf27 
$ ipython
(Python 2.7)
> run mk_bias.py
```

After running this code, please check if ``Mbias.fits`` is created well in `./calibrations/`.
* Running time: ~240 sec (for 15 raw bias files), ~800 sec (for 30 files)
