#!/usr/bin/env python3

from pathlib import Path
import time

import numpy as np
import glob, os, copy
from astropy.io import fits
from astropy.stats import sigma_clip
import shutil
ds9_path = shutil.which('ds9')
import subprocess

import init_cfg as ic


# ----- Functions ----- #
def wavelength_axis(header):
    """Return the linear spectral-coordinate array."""

    n_wave = header["NAXIS3"]
    pixel = np.arange(n_wave, dtype=np.float64) + 1.0

    if "CD3_3" in header:
        dw = header["CD3_3"]
    else:
        dw = header["CDELT3"]

    wave = header["CRVAL3"] + (pixel - header["CRPIX3"]) * dw
    return wave


def spectral_regrid_linear(wave, sci, var, wave_new):
    """
    Linear interpolation along axis=0.
    The variance is propagated as
        var_new = w0**2 * var0 + w1**2 * var1
    assuming adjacent input wavelength bins are independent.
    This function assumes SCI represents a flux-density-like quantity
    sampled at wavelength-bin centers.
    """

    wave = np.asarray(wave, dtype=np.float64)

    # Ensure increasing wavelength
    if not np.array_equal(wave, np.sort(wave)):
        wave = wave[::-1]
        sci = sci[::-1]
        var = var[::-1]

    assert np.all(np.diff(wave) > 0)
    
    # Interpolation index: wave[lo] <= wave_new <= wave[hi]
    hi = np.searchsorted(wave, wave_new, side="left")

    hi = np.clip(hi, 1, wave.size - 1)
    lo = hi - 1

    denominator = wave[hi] - wave[lo]
    frac = (wave_new - wave[lo]) / denominator    # Grid fraction within a wavelength bin

    a0 = (1.0 - frac)[:, None, None]
    a1 = frac[:, None, None]

    sci0 = np.asarray(sci[lo], dtype=np.float64)    # Lower-index science data
    sci1 = np.asarray(sci[hi], dtype=np.float64)    # Upper-index science data

    var0 = np.asarray(var[lo], dtype=np.float64)    # Lower-index variance data
    var1 = np.asarray(var[hi], dtype=np.float64)    # Upper-index variance data

    good0 = (np.isfinite(sci0) & np.isfinite(var0) & (var0 >= 0))
    good1 = (np.isfinite(sci1) & np.isfinite(var1) & (var1 >= 0))
    
    # Weighted interpolation for new wavelength bins
    # If one of the adjacent bins is invalid, renormalize using the
    # remaining valid bin.
    weight0 = a0 * good0
    weight1 = a1 * good1
    weight_sum = weight0 + weight1

    sci_numerator = (np.where(good0, a0 * sci0, 0.0) + np.where(good1, a1 * sci1, 0.0))
    var_numerator = (np.where(good0, a0**2 * var0, 0.0) + np.where(good1, a1**2 * var1, 0.0))

    sci_new = np.full((wave_new.size, sci.shape[1], sci.shape[2]),
                      np.nan, dtype=np.float64)
    var_new = np.full_like(sci_new, np.nan)

    np.divide(sci_numerator, weight_sum,
              out=sci_new,
              where=weight_sum > 0)

    np.divide(var_numerator, weight_sum**2,
              out=var_new,
              where=weight_sum > 0)

    inside = ((wave_new >= wave[0]) & \
              (wave_new <= wave[-1]))

    sci_new[~inside] = np.nan
    var_new[~inside] = np.nan

    return sci_new.astype(np.float32), var_new.astype(np.float32)


def spatial_shift_bilinear(sci, var, shift_y, shift_x):
    """
    Apply a constant bilinear spatial shift to a full 3D cube.
    The shift convention follows scipy.ndimage.shift:
        positive shift_y -> source moves toward larger Y indices
        positive shift_x -> source moves toward larger X indices
    Variance is propagated using squared interpolation coefficients.
    """

    n_wave, ny, nx = sci.shape

    # Input coordinate sampled by each output coordinate.
    source_y = np.arange(ny, dtype=np.float64) - shift_y
    source_x = np.arange(nx, dtype=np.float64) - shift_x

    y0 = np.floor(source_y).astype(int)
    x0 = np.floor(source_x).astype(int)

    y1 = y0 + 1
    x1 = x0 + 1

    fy = source_y - y0
    fx = source_x - x0

    ### corners = [index_y, weight_y, index_x, weight_x]
    ###  with (index_y, index_x) are the target-frame coordinates.
    corners = (
        (y0, 1.0 - fy, x0, 1.0 - fx),
        (y0, 1.0 - fy, x1, fx),
        (y1, fy,       x0, 1.0 - fx),
        (y1, fy,       x1, fx),
    )

    sci_numerator = np.zeros(sci.shape, dtype=np.float64)
    var_numerator = np.zeros(sci.shape, dtype=np.float64)
    weight_sum    = np.zeros(sci.shape, dtype=np.float64)

    for iy, wy, ix, wx in corners:
        coordinate_valid = (
            (iy[:, None] >= 0) & (iy[:, None] < ny) & \
            (ix[None, :] >= 0) & (ix[None, :] < nx)
        )

        iy_safe = np.clip(iy, 0, ny - 1)
        ix_safe = np.clip(ix, 0, nx - 1)

        corner_sci = sci[:, iy_safe[:, None], ix_safe[None, :]]
        corner_var = var[:, iy_safe[:, None], ix_safe[None, :]]
        spatial_weight = (wy[:, None] * wx[None, :])

        good = (coordinate_valid[None, :, :] & \
                np.isfinite(corner_sci) & \
                np.isfinite(corner_var) & \
                (corner_var >= 0))

        effective_weight = (spatial_weight[None, :, :] * good)

        sci_numerator += (effective_weight * np.where(good, corner_sci, 0.0))
        var_numerator += (effective_weight**2 * np.where(good, corner_var, 0.0))
        weight_sum    += effective_weight

    sci_shifted = np.full(sci.shape, np.nan, dtype=np.float64)
    var_shifted = np.full(sci.shape, np.nan, dtype=np.float64)

    np.divide(sci_numerator, weight_sum, out=sci_shifted, where=weight_sum > 0)
    np.divide(var_numerator, weight_sum**2, out=var_shifted, where=weight_sum > 0)

    # For fully covered interior pixels this is approximately 1.
    coverage = weight_sum.astype(np.float32)

    return (
        sci_shifted.astype(np.float32),
        var_shifted.astype(np.float32),
        coverage,
    )


def combine_weighted(science_stack, variance_stack,
                     sigma=3.0, maxiters=5):
    """
    Sigma-clipped inverse-variance weighted combination.
    Inputs have shape
        (n_cube, n_wave, ny, nx)
    """
    
    base_mask = (~np.isfinite(science_stack) | \
                 ~np.isfinite(variance_stack) | \
                 (variance_stack <= 0))

    # Sigma clipping is not particularly robust when only a few exposures are available.
    if science_stack.shape[0] >= 5 and sigma is not None:
        clipped = sigma_clip(
            np.ma.array(science_stack, mask=base_mask),
            sigma=sigma,
            maxiters=maxiters,
            cenfunc="median",
            stdfunc="mad_std",
            axis=0,
            masked=True,
        )

        final_mask = np.ma.getmaskarray(clipped)    # 'True' means masked values.
    else:
        final_mask = base_mask

    weights = np.zeros(variance_stack.shape, dtype=np.float64)
    np.divide(1.0, variance_stack, out=weights, where=~final_mask)

    science_valid = np.where(final_mask, 0.0, science_stack)

    numerator = np.sum(weights * science_valid, axis=0, dtype=np.float64)
    weight_sum = np.sum(weights, axis=0, dtype=np.float64)

    combined_science = np.full(science_stack.shape[1:], np.nan, dtype=np.float64)
    combined_variance = np.full_like(combined_science, np.nan)

    np.divide(
        numerator,
        weight_sum,
        out=combined_science,
        where=weight_sum > 0,
    )

    np.divide(
        1.0,
        weight_sum,
        out=combined_variance,
        where=weight_sum > 0,
    )

    n_contributing = np.sum(
        ~final_mask,
        axis=0,
    ).astype(np.uint16)

    return (
        combined_science.astype(np.float32),
        combined_variance.astype(np.float32),
        n_contributing,
    )


def update_spectral_header(header, wave_new):
    header = header.copy()

    header["CRPIX3"] = 1.0
    header["CRVAL3"] = float(wave_new[0])

    dw = float(wave_new[1] - wave_new[0])

    if "CD3_3" in header:
        header["CD3_3"] = dw

    if "CDELT3" in header or "CD3_3" not in header:
        header["CDELT3"] = dw

    return header


# ----- Main workflows ----- #
# def main():
start_time = time.time()

cube_paths = [Path(name) for name in ic.cube_list]
# cube_paths = cube_paths[0:1]
# cube_paths  = [Path(name) for name in sorted(glob.glob("cstxeqxbrg*_3D.fits"))]

offset_file = "offset.txt"
offset_x, offset_y, offset_dist = np.loadtxt(offset_file, unpack=True)
# offset_x, offset_y, offset_dist = offset_x[0:1], offset_y[0:1], offset_dist[0:1]

assert len(cube_paths) == len(offset_x)

### Determine the common wavelength range for all cubes
wavelength_ranges = []

for path in cube_paths:
    with fits.open(path, memmap=True) as hdul:
        wave = wavelength_axis(hdul[1].header)
    wavelength_ranges.append((np.nanmin(wave), np.nanmax(wave)))

common_start = max(item[0] for item in wavelength_ranges)
common_end   = min(item[1] for item in wavelength_ranges)

margin = 10.0    # Optional edge margin [AA]
wave_start = common_start + margin
wave_end   = common_end   - margin

assert wave_start < wave_end

wave_new = np.arange(wave_start, wave_end+0.5*ic.wav_intv, ic.wav_intv,
                     dtype=np.float64)


### Initial mask of each cube
aligned_science  = []
aligned_variance = []

primary_header = None
reference_sci_header = None
reference_var_header = None

for index, (path, dx, dy) in enumerate(
    zip(cube_paths, offset_x, offset_y)):
    print(f"[{index + 1}/{len(cube_paths)}] {path.name}")

    with fits.open(path, memmap=True) as hdul:
        if index == 0:
            primary_header = hdul[0].header.copy()
            reference_sci_header = hdul[1].header.copy()
            reference_var_header = hdul[2].header.copy()

        sci  = np.array(hdul[1].data, dtype=np.float32, copy=True)
        var  = np.array(hdul[2].data, dtype=np.float32, copy=True)
        wave = wavelength_axis(hdul[1].header)

    # A small guard region needed for interpolation
    native_dw = np.nanmedian(np.abs(np.diff(wave)))
    guard = 2.0 * native_dw

    use = ((wave >= wave_start - guard) & \
           (wave <= wave_end + guard))

    assert np.sum(use) >= 2

    sci, var, wave = sci[use], var[use], wave[use]

    # Spectral alignment
    sci, var = spectral_regrid_linear(wave, sci, var, wave_new)

    # Match the original code's sign convention:
    # shift=(-offset_Y, -offset_X)
    sci, var, coverage = spatial_shift_bilinear(
        sci, var,
        shift_y=-dy,
        shift_x=-dx,
    )

    # Reject partially covered boundary pixels.
    good_coverage = coverage >= 0.999

    sci[~good_coverage] = np.nan
    var[~good_coverage] = np.nan
    
    aligned_science.append(sci)
    aligned_variance.append(var)    

science_stack  = np.stack(aligned_science, axis=0)
variance_stack = np.stack(aligned_variance, axis=0)

'''
##### Check #####
sci_header = update_spectral_header(reference_sci_header, wave_new)
var_header = update_spectral_header(reference_var_header, wave_new)

for index in range(len(cube_paths)):
    sci_hdu  = fits.ImageHDU(data=science_stack[index, :, :, :], header=sci_header, name="SCI")
    var_hdu  = fits.ImageHDU(data=variance_stack[index, :, :, :], header=var_header, name="VAR")

    output_path = Path(ic.dir_cmb) / f"a{str(cube_paths[index]).split("/")[-1].split('cstxeqxbrg')[-1]}"

    new_hdul = fits.HDUList([fits.PrimaryHDU(header=primary_header), sci_hdu, var_hdu])
    new_hdul.writeto(output_path, overwrite=True)
#################
'''

final_science, final_variance, n_contributing = (
    combine_weighted(
        science_stack,
        variance_stack,
        sigma=None,
        maxiters=None,
    )
)

sci_header = update_spectral_header(reference_sci_header, wave_new)
var_header = update_spectral_header(reference_var_header, wave_new)


# The aligned product is on the reference cube's spatial WCS.
sci_hdu  = fits.ImageHDU(data=final_science , header=sci_header, name="SCI")
var_hdu  = fits.ImageHDU(data=final_variance, header=var_header, name="VAR")
nexp_hdu = fits.ImageHDU(data=n_contributing, header=sci_header, name="NEXP")

output_path = Path(ic.dir_cmb) / "combined_cube_3D.fits"

new_hdul = fits.HDUList([fits.PrimaryHDU(header=primary_header), sci_hdu, var_hdu, nexp_hdu])
new_hdul.writeto(output_path, overwrite=True)

print(f"Written: {output_path}")
print(
    f"Wavelength range: "
    f"{wave_new[0]:.2f} - {wave_new[-1]:.2f} AA"
)
print(f"Number of wavelength bins: {wave_new.size}")
print(f"Elapsed time: {time.time() - start_time:.2f} s")



# if __name__ == "__main__":
    # main()


