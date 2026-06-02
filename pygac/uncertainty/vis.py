#!/usr/bin/env python

# Copyright (c) 2014-2015, 2019 Pytroll Developers

# Author(s):

#   Nicole Yaghnam <nicole.yaghnam@npl.co.uk>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Uncertainty information based on the noaa.py calibration for VIS channels
"""
from __future__ import division

import numpy as np
import xarray as xr

from pygac.calibration.noaa import Calibrator
from pygac.uncertainty.ir import allan_deviation, get_bad_space_counts, get_uncert_parameter_thresholds


def _vis_channel_noise(space_2d, bad_scans, window):
    """Compute noise and averaged-noise scalars for one VIS channel.

    Parameters
    ----------
    space_2d : (N, P) array
        Space-view counts for this channel.
    bad_scans : (N,) int8 array
        Bad-scan mask from :func:`_compute_vis_bad_scans`.
    window : int
        Averaging window size (number of scanlines in the kernel).

    Returns
    -------
    noise : scalar
        Allan deviation combined with digitisation uncertainty.
    av_noise : scalar
        Noise after averaging over *window* × 10 measurements.
    """
    raw = allan_deviation(space_2d, bad_scan=bad_scans)
    noise = np.sqrt(raw ** 2 + 1.0 / 3)
    av_noise = noise / np.sqrt(window * 10)
    return noise, av_noise


def _compute_vis_bad_scans(total_space, chan_3a):
    """Flag scanlines that contain bad space-view counts in any VIS channel.

    Parameters
    ----------
    total_space : (N, P, C) array
        Space-view counts; C >= 2 (C=3 when *chan_3a* is True).
    chan_3a : bool
        Whether the 1.6 µm (channel 3a) space-view data is present.

    Returns
    -------
    bad_scans : (N,) int8 array
        1 for each scanline with at least one bad pixel in any active channel.
    """
    n_chans = 3 if chan_3a else 2
    bad_per_pixel = np.zeros(total_space.shape[:2], dtype=bool)
    for k in range(n_chans):
        bad_per_pixel |= get_bad_space_counts(total_space[:, :, k])
    return bad_per_pixel.any(axis=1).astype(np.int8)


def get_noise(total_space, window, chan_3a):
    """Get noise estimates from the counts.

    .. deprecated::
        Use :func:`_compute_vis_bad_scans` and :func:`_vis_channel_noise` instead.
    """
    import warnings
    warnings.warn(
        "get_noise is deprecated; use _compute_vis_bad_scans + _vis_channel_noise instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    bad_scans = _compute_vis_bad_scans(total_space, chan_3a)
    noise1, av_noise1 = _vis_channel_noise(total_space[:, :, 0], bad_scans, window)
    noise2, av_noise2 = _vis_channel_noise(total_space[:, :, 1], bad_scans, window)
    if chan_3a:
        noise3, av_noise3 = _vis_channel_noise(total_space[:, :, 2], bad_scans, window)
    else:
        noise3, av_noise3 = None, None
    return noise1, noise2, noise3, av_noise1, av_noise2, av_noise3, bad_scans


def _vis_random_uncert(noise, av_noise, gain, counts, mean_space):
    """Random radiance uncertainty for one VIS channel scanline.

    Parameters
    ----------
    noise : scalar
        Space-view noise (Allan dev + digitisation).
    av_noise : scalar
        Noise after averaging over the kernel window.
    gain : scalar
        Calibration slope.
    counts : (P,) array
        Per-pixel counts.
    mean_space : scalar
        Mean space-view count for this scanline.

    Returns
    -------
    uncert : (P,) array
        Random radiance uncertainty.
    Rcal : (P,) array
        Calibrated radiance.
    """
    Rcal = gain * (counts - mean_space)
    uncert = np.sqrt((gain * av_noise) ** 2 + (gain * noise) ** 2)
    return np.full(counts.shape, uncert), Rcal


def get_random(noise, av_noise, gain, cal, year, jday, C, D):
    """Get the random parts of the vis calibration uncertainty.

    Done per scanline"""
    #
    # Get time since launch in years
    #
    # l_date = Calibrator.date2float(cal.date_of_launch)
    # t = (year + jday/365.0) - l_date
    #
    # Measurement Function
    #
    Rcal = gain*(C-D)
    #
    # Gain part for all noise sources
    #
    dRcal_dC = gain
    dRcal_dD = -gain

    uncert = (dRcal_dD**2)*(av_noise**2) + (dRcal_dC**2)*(noise**2)

    return np.sqrt(uncert), Rcal

def get_reflectance(Rcal, d_se, sza):
    refl = (Rcal * d_se ** 2) / np.cos(sza) / 10

    return refl

_U_SYS_BASE = np.sqrt(0.025**2 + 0.025**2 + 0.01**2 + 0.015**2 + 0.02**2 + 0.025**2)
_U_WV = 0.015


def _vis_sys_uncert(counts, mean_space, gain, include_water_vapour):
    """Systematic radiance uncertainty for one VIS channel.

    Parameters
    ----------
    counts : (P,) array
        Per-pixel counts.
    mean_space : scalar
        Mean space-view count for this scanline.
    gain : scalar
        Calibration slope.
    include_water_vapour : bool
        True for the 0.86 µm channel, which has additional WV uncertainty.

    Returns
    -------
    uncert : (P,) array
        Systematic radiance uncertainty [counts × gain].
    """
    usys = np.sqrt(_U_SYS_BASE**2 + _U_WV**2) if include_water_vapour else _U_SYS_BASE
    usys_total = usys * gain
    dRcal_dS = counts - mean_space
    return np.sqrt(dRcal_dS**2 * usys_total**2)


def get_sys(channel, C, D, gain):
    """Get the systematic parts of the vis calibration uncertainty.

    Comprised of:
        MODIS Reflectance uncertainty - 2.5%
        SBAF Correction uncertainty (MODIS SNOs) - 2.5%
        Dome-C Surface Reflectance variation - 1%
        Temporal Stability - 1.5%
        Libya Surface Reflectance variation - 2%
        SBAF Correction uncertainty (PICS) - 2.5%
        Water Vapour effects (Channel 2 only) - 1.5%
    """
    dRcal_dS = (C-D)
    usys = np.sqrt(0.025**2 + 0.025**2 + 0.01**2 + 0.015**2 + 0.02**2 + 0.025**2)
    # If channel = 2, add water vapour uncertainty
    #
    U_WV = 0.015
    if channel == 2:
        usys_tot = np.sqrt(usys**2 + U_WV**2)
    else:
        usys_tot = usys

    usys_tot *= gain

    uncert = (dRcal_dS**2)*(usys_tot**2)

    return np.sqrt(uncert)

def get_vars(ds, channel):
    """Get variables from xarray"""

    mean_space = ds["full_space_counts"].values[:, :, channel].mean(axis=1)
    counts = ds["counts"].values[:, :, channel]

    return mean_space, counts

def get_gain(s0, s1, s2, t, cal, channel):
    if np.isnan(cal.gain_switch).all():
        glow = ghigh = np.ones(3)
    else:
        glow = np.array([0.5, 0.5, 0.25])
        ghigh = np.array([1.5, 1.5, 1.75])

    s0_l = s0 * glow[channel]
    s0_h = s0 * ghigh[channel]

    stl = s0_l*(100 + s1*t + s2*t**2)/100
    sth = s0_h*(100 + s1*t + s2*t**2)/100

    gain = (stl + sth)/2

    return gain


def vis_uncertainty(ds, mask):
    """Create the uncertainty components for the vis channels. These include

    1) Random
        a) Noise
        b) Digitisation
    2) Systematic
        a) MODIS Reflectance uncertainty
        b) SBAF Correction uncertainty
        c) Water Vapour effects
        d) Surface Reflectance variation
        e) Temporal Stability

    Inputs:
          ds : Input xarray dataset containing data for calibration
        mask : pygac mask from reader
    Outputs:
      uncert : xarray dataset containing random and systematic uncertainty
               components
    """
    #
    # Define averaging kernel based on value in noaa.py
    # Also get the thresholds for solar contamination detection
    #
    window, solar_contam_threshold, sza_threshold = \
        get_uncert_parameter_thresholds(vischans=True)

    # if ds["channels"].values.shape[1] == 409:
    #     gacdata = True
    # else:
    #     gacdata = False

    # avhrr_name = ds.attrs["spacecraft_name"]

    #
    # Test for channel 3a
    # Testing for the 3a/3b toggle needs to be added at the scanline level
    #
    #if KLMReader._get_vis_channels_to_calibrate == [0,1,2]:
    if ds["channels"].values.shape[2] == 6:
        if np.all(ds["channels"].values[:,:,2] == np.nan):
            chan_3a = False
        else:
            chan_3a = True
    else:
        chan_3a = False

    #
    # Get calibration coefficients
    #
    cal = Calibrator(
        ds.attrs["spacecraft_name"])
    s0_1 = cal.s0[0]
    s0_2 = cal.s0[1]
    s0_3 = cal.s0[2]
    s1_1 = cal.s1[0]
    s1_2 = cal.s1[1]
    s1_3 = cal.s1[2]
    s2_1 = cal.s2[0]
    s2_2 = cal.s2[1]
    s2_3 = cal.s2[2]

    #
    # Get time since launch in years
    #
    times = ds.coords["times"]
    start_time = times[0].dt
    year = start_time.year.item()
    jday = start_time.dayofyear.item()

    #
    # Get variables for 10 sampled case
    #
    nb_refl = 3 if chan_3a else 2
    total_space = ds["full_space_counts"].values[:, :, nb_refl:]

    #
    # Noise elements
    #
    bad_scan = _compute_vis_bad_scans(total_space, chan_3a)
    noise1, av_noise1 = _vis_channel_noise(total_space[:, :, 0], bad_scan, window)
    noise2, av_noise2 = _vis_channel_noise(total_space[:, :, 1], bad_scan, window)
    if chan_3a:
        noise3, av_noise3 = _vis_channel_noise(total_space[:, :, 2], bad_scan, window)

    #
    # Get variables used on the calibration
    #
    mean_space_1, counts_1 = get_vars(ds, 0)
    mean_space_2, counts_2 = get_vars(ds, 1)
    if chan_3a:
        mean_space_3a, counts_3a = get_vars(ds, 2)

    #
    # Calibration slope
    #
    l_date = Calibrator.date2float(cal.date_of_launch)
    t = (year + jday / 365.0) - l_date
    gain_1 = get_gain(s0_1, s1_1, s2_1, t, cal, 0)
    gain_2 = get_gain(s0_2, s1_2, s2_2, t, cal, 1)
    if chan_3a:
        gain_3 = get_gain(s0_3, s1_3, s2_3, t, cal, 2)

    #
    # Vectorised random uncertainty — constant per channel across all scanlines
    #
    rand_const_63 = np.sqrt((gain_1 * av_noise1) ** 2 + (gain_1 * noise1) ** 2)
    rand_const_86 = np.sqrt((gain_2 * av_noise2) ** 2 + (gain_2 * noise2) ** 2)

    rcal_rand_63 = np.full(counts_1.shape, rand_const_63)
    rcal_rand_86 = np.full(counts_2.shape, rand_const_86)
    if chan_3a:
        rand_const_12 = np.sqrt((gain_3 * av_noise3) ** 2 + (gain_3 * noise3) ** 2)
        rcal_rand_12 = np.full(counts_3a.shape, rand_const_12)
    else:
        rcal_rand_12 = np.full(counts_2.shape, np.nan)

    #
    # Vectorised systematic uncertainty — varies per pixel (counts - mean_space)
    #
    rcal_sys_63 = _vis_sys_uncert(counts_1, mean_space_1[:, np.newaxis], gain_1, include_water_vapour=False)
    rcal_sys_86 = _vis_sys_uncert(counts_2, mean_space_2[:, np.newaxis], gain_2, include_water_vapour=True)
    if chan_3a:
        rcal_sys_12 = _vis_sys_uncert(counts_3a, mean_space_3a[:, np.newaxis], gain_3, include_water_vapour=False)
    else:
        rcal_sys_12 = np.full(counts_2.shape, np.nan)

    #
    # Rcal for solar contamination detection (bad scans → 0.0 to match legacy)
    #
    Rcal_1 = gain_1 * (counts_1 - mean_space_1[:, np.newaxis])
    Rcal_1[bad_scan == 1] = 0.0

    #
    # Apply bad-scan mask
    #
    bad = bad_scan == 1
    for arr in (rcal_rand_63, rcal_rand_86, rcal_rand_12,
                rcal_sys_63, rcal_sys_86, rcal_sys_12):
        arr[bad] = np.nan

    # Solar contamination check
    d_se = ds.attrs["sun_earth_distance_correction_factor"]
    sza = ds["sun_zen"].values
    refl_1 = get_reflectance(Rcal_1, d_se, sza)
    contam_pixels = np.zeros(refl_1.shape, dtype=np.int8)
    gd = (refl_1 > solar_contam_threshold) & (sza > sza_threshold)
    if np.sum(gd) > 0:
        contam_pixels[gd] = 1

    #
    # Assemble output dataset
    #
    random = np.stack([rcal_rand_63, rcal_rand_86, rcal_rand_12], axis=-1)
    systematic = np.stack([rcal_sys_63, rcal_sys_86, rcal_sys_12], axis=-1)

    time = (ds["times"].values - np.datetime64("1970-01-01 00:00:00")) / np.timedelta64(1, "s")

    return xr.Dataset({
        "times": xr.DataArray(time, dims=["times"],
                              attrs={"long_name": "scanline time", "units": "seconds since 1970-01-01"}),
        "across_track": xr.DataArray(np.arange(random.shape[1]), dims=["across_track"]),
        "vis_channels": xr.DataArray(np.array([1, 2, 3]), dims=["vis_channels"]),
        "random": xr.DataArray(random, dims=["times", "across_track", "vis_channels"],
                               attrs={"long_name": "Random uncertainties", "units": ""}),
        "systematic": xr.DataArray(systematic, dims=["times", "across_track", "vis_channels"],
                                   attrs={"long_name": "Systematic uncertainties", "units": ""}),
        "solar_fov_contam": xr.DataArray(
            contam_pixels, dims=["times", "across_track"],
            attrs={"long_name": "Flag for in FOV solar contamination (0=none, 1=contaminated)",
                   "units": ""}),
    })
