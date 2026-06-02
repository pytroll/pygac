#!/usr/bin/env python

# Copyright (c) 2014-2015, 2019 Pytroll Developers

# Author(s):

#   Jonathan Mittaz <j.mittaz@reading.ac.uk>

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

"""Uncertainty information based on the noaa.py calibration for IR channels
"""
from __future__ import division

from contextlib import contextmanager
from dataclasses import dataclass
from tempfile import gettempdir

import numpy as np
import xarray as xr
from scipy.optimize import curve_fit

from pygac.calibration.noaa import Calibrator, get_prt_nos

#: Platforms that only have 3.7 µm and 11 µm IR channels (no 12 µm).
_NO_TWELVE_MICRON = frozenset({"tirosn", "noaa06", "noaa08", "noaa10"})


@dataclass(frozen=True)
class IRChannelSpec:
    """Immutable identity and calibration constants for one IR channel.

    Three independent index systems exist in the IR calibration code:

    * ``cal_index``     — 0/1/2, indexes into :class:`~pygac.calibration.noaa.Calibrator` arrays.
    * ``output_index``  — 3/4/5, the channel slot in the merged xr.Dataset output.
    * ``helper_channel``— 1/2/3, the legacy positional argument expected by
                          :func:`get_random` and :func:`get_sys`.
    """

    cal_index: int
    output_index: int
    helper_channel: int
    label: str
    conv: object          # convBT — not hashable, so frozen-dataclass equality is identity-based
    space_radiance: float
    nonlin_coeffs: tuple  # (c0, c1, c2); (0., 0., 0.) for 3.7 µm
    has_nonlinear: bool   # False for 3.7 µm only


def ir_channel_specs(cal, platform):
    """Return the :class:`IRChannelSpec` list for *platform*.

    Parameters
    ----------
    cal : Calibrator
        Calibrator instance for the platform.
    platform : str
        Spacecraft name, e.g. ``"noaa14"``.

    Returns
    -------
    list[IRChannelSpec]
        Two specs for tirosn/noaa06/noaa08/noaa10 (no 12 µm channel),
        three specs for all other platforms.
    """
    specs = [
        IRChannelSpec(
            cal_index=0,
            output_index=3,
            helper_channel=1,
            label="3.7",
            conv=convBT(cal, 0),
            space_radiance=0.0,
            nonlin_coeffs=(0.0, 0.0, 0.0),
            has_nonlinear=False,
        ),
        IRChannelSpec(
            cal_index=1,
            output_index=4,
            helper_channel=2,
            label="11",
            conv=convBT(cal, 1),
            space_radiance=cal.space_radiance[1],
            nonlin_coeffs=(cal.b[1, 0], cal.b[1, 1], cal.b[1, 2]),
            has_nonlinear=True,
        ),
    ]
    if platform not in _NO_TWELVE_MICRON:
        specs.append(
            IRChannelSpec(
                cal_index=2,
                output_index=5,
                helper_channel=3,
                label="12",
                conv=convBT(cal, 2),
                space_radiance=cal.space_radiance[2],
                nonlin_coeffs=(cal.b[2, 0], cal.b[2, 1], cal.b[2, 2]),
                has_nonlinear=True,
            )
        )
    return specs


@dataclass
class IRTelemetry:
    """Channel-independent telemetry computed once for the whole orbit.

    The PRT pipeline (mean_prt_counts → temperatures) is identical across
    all three IR channels. Computing it once avoids the 3× redundancy that
    was present in the original code.
    """

    tict: np.ndarray   # smoothed ICT temperature (scanlines,)
    prt1: np.ndarray   # per-PRT smoothed temperatures (scanlines,)
    prt2: np.ndarray
    prt3: np.ndarray
    prt4: np.ndarray


def extract_ir_telemetry(ds, cal, mask, window, prt_threshold, gacdata):
    """Compute the channel-independent PRT telemetry for the orbit.

    This is the channel-independent subset of the former :func:`get_vars`
    logic: it reads ``mean_prt_counts``, interpolates bad PRT readings,
    converts counts → temperatures, and applies the smoothing convolution.

    Parameters
    ----------
    ds : xr.Dataset
    cal : Calibrator
    mask : np.ndarray  (boolean, shape=(scanlines,))
    window : int
    prt_threshold : float
    gacdata : bool

    Returns
    -------
    IRTelemetry
    """
    line_numbers = ds["scan_line_index"].data
    prt = ds["mean_prt_counts"].values.copy()

    gd = ~np.isfinite(prt)
    if np.sum(gd) > 0:
        prt[gd] = 0

    iprt = get_prt_nos(prt, prt_threshold, line_numbers, gacdata)

    for prt_idx in range(1, 5):
        ifix = np.where(np.logical_and(iprt == prt_idx, prt <= prt_threshold))
        if len(ifix[0]):
            inofix = np.where(np.logical_and(iprt == prt_idx, prt > prt_threshold))
            if len(inofix[0]):
                prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])
            else:
                raise IndexError(f"No good prt{prt_idx} data")

    tprt = np.polynomial.polynomial.polyval(prt, cal.d[:, iprt], tensor=False)

    weighting_function = np.ones(window, dtype=float) / window
    half = (window - 1) // 2
    half_up = (window + 1) // 2

    def _convolve_clamp(arr):
        out = np.convolve(arr, weighting_function, "same")
        out[:half] = out[half]
        out[-half:] = out[-half_up]
        return out

    tprt_interp = np.copy(tprt)
    zeros = iprt == 0
    nonzeros = ~zeros
    tprt_interp[zeros] = np.interp(zeros.nonzero()[0], nonzeros.nonzero()[0], tprt[nonzeros])
    tict = _convolve_clamp(tprt_interp)

    def _per_prt(idx):
        arr = np.copy(tprt)
        mask_ = (iprt == 0) | (iprt != idx)
        nonmask = ~mask_
        arr[mask_] = np.interp(mask_.nonzero()[0], nonmask.nonzero()[0], tprt[nonmask])
        return _convolve_clamp(arr)

    return IRTelemetry(
        tict=tict,
        prt1=_per_prt(1),
        prt2=_per_prt(2),
        prt3=_per_prt(3),
        prt4=_per_prt(4),
    )


@dataclass
class IRChannelData:
    """Per-channel arrays produced by the IR calibration pre-processing step.

    All arrays are indexed by scanline (1-D) except ``ce`` which is
    (scanlines × pixels) 2-D. These are mutable: the per-scanline loop in
    the orchestrator reads from them in place.
    """

    spec: IRChannelSpec
    noise: np.ndarray
    av_noise: np.ndarray
    av_ict_noise: np.ndarray
    cs: np.ndarray
    cict: np.ndarray
    ce: np.ndarray


def build_ir_channel_data(
    ds,
    specs,
    total_space,
    total_ict,
    window,
    prt_threshold,
    ict_threshold,
    space_threshold,
    gacdata,
    cal,
    mask,
):
    """Build one :class:`IRChannelData` per channel spec and return shared telemetry.

    Parameters
    ----------
    ds : xr.Dataset
        Calibrated input dataset.
    specs : list[IRChannelSpec]
        Channel specs from :func:`ir_channel_specs`.
    total_space, total_ict : np.ndarray
        Raw space/ICT count arrays ``(scanlines, samples_per_line, channels)``.
    window, prt_threshold, ict_threshold, space_threshold : float
        Uncertainty parameter thresholds from :func:`get_uncert_parameter_thresholds`.
    gacdata : bool
        True for GAC (409-pixel) data.
    cal : Calibrator
        Calibrator instance for the platform.
    mask : np.ndarray
        Boolean scanline mask.

    Returns
    -------
    channels : list[IRChannelData]
        One entry per spec, in spec order.
    bad_scan : np.ndarray
        1-D bool-like array flagging bad scanlines (from :func:`get_noise`).
    Tict : np.ndarray
        Smoothed ICT temperature (shared across channels; same PRT pipeline).
    ict1, ict2, ict3, ict4 : np.ndarray
        Per-PRT smoothed temperatures (needed by :func:`get_gainval`).
    """
    noise_all, bad_scan = _compute_noise_arrays(specs, total_space, total_ict, window)

    # PRT pipeline runs only once — it is channel-independent.
    telemetry = extract_ir_telemetry(ds, cal, mask, window, prt_threshold, gacdata)
    Tict = telemetry.tict
    ict1, ict2, ict3, ict4 = telemetry.prt1, telemetry.prt2, telemetry.prt3, telemetry.prt4

    channels = []
    for spec in specs:
        cs, cict, ce = _get_channel_arrays(
            ds, spec.cal_index, Tict, mask, ict_threshold, space_threshold, window,
        )
        noise, av_noise, av_ict_noise = noise_all[spec.cal_index]
        channels.append(IRChannelData(
            spec=spec,
            noise=noise,
            av_noise=av_noise,
            av_ict_noise=av_ict_noise,
            cs=cs,
            cict=cict,
            ce=ce,
        ))

    return channels, bad_scan, Tict, ict1, ict2, ict3, ict4


def _channel_noise(space_2d, ict_2d, bad_scans, window):
    """Compute noise estimates for a single IR channel.

    Parameters
    ----------
    space_2d : np.ndarray, shape (scanlines, counts_per_line)
    ict_2d   : np.ndarray, shape (scanlines, counts_per_line)
    bad_scans : np.ndarray, shape (scanlines,), dtype int8
    window : int  Number of scanlines in the smoothing window.

    Returns
    -------
    noise : float  Space-count noise (Allan deviation + digitisation).
    av_noise : float  Noise after averaging over window * 10 measurements.
    av_ict_noise : float  ICT noise after averaging.
    """
    noise = np.sqrt(allan_deviation(space_2d, bad_scan=bad_scans) ** 2 + 1.0 / 3)
    ict_noise = np.sqrt(allan_deviation(ict_2d, bad_scan=bad_scans) ** 2 + 1.0 / 3)
    sqrt_window = np.sqrt(window * 10)
    return noise, noise / sqrt_window, ict_noise / sqrt_window


def _compute_bad_scans(specs, total_space, total_ict):
    """Flag scanlines that have bad space-count data in any channel.

    Parameters
    ----------
    specs : list[IRChannelSpec]
    total_space : np.ndarray, shape (scanlines, counts_per_line, 3)
    total_ict   : np.ndarray, shape (scanlines, counts_per_line, 3)

    Returns
    -------
    bad_scans : np.ndarray, shape (scanlines,), dtype int8
        1 where any channel reports bad data, 0 elsewhere.
    """
    bad_per_channel = [
        get_bad_space_counts(
            total_space[:, :, s.cal_index],
            ict_data=total_ict[:, :, s.cal_index],
        )
        for s in specs
    ]
    stacked = np.stack(bad_per_channel, axis=0)  # (n_channels, scanlines, counts)
    return stacked.any(axis=(0, 2)).astype(np.int8)


def _compute_noise_arrays(specs, total_space, total_ict, window):
    """Compute per-channel noise estimates.

    Parameters
    ----------
    specs : list[IRChannelSpec]
    total_space : np.ndarray, shape (scanlines, counts_per_line, 3)
    total_ict   : np.ndarray, shape (scanlines, counts_per_line, 3)
    window : int

    Returns
    -------
    noise_by_cal_index : dict[int, tuple[float, float, float]]
        Maps cal_index → (noise, av_noise, av_ict_noise).
    bad_scans : np.ndarray, shape (scanlines,), dtype int8
    """
    bad_scans = _compute_bad_scans(specs, total_space, total_ict)
    return {
        s.cal_index: _channel_noise(
            total_space[:, :, s.cal_index],
            total_ict[:, :, s.cal_index],
            bad_scans,
            window,
        )
        for s in specs
    }, bad_scans


def _get_channel_arrays(ds, cal_index, tict, mask, ict_threshold, space_threshold, window):
    """Extract the channel-specific calibration arrays from *ds*.

    This is the channel-specific subset of the former :func:`get_vars` logic:
    it reads space/ICT/earth counts for *cal_index*, interpolates over bad
    values, and applies the smoothing convolution.  The PRT-derived *tict*
    (smoothed ICT temperature) is passed in from :func:`extract_ir_telemetry`
    so it is not recomputed per channel.

    Parameters
    ----------
    ds : xr.Dataset
    cal_index : int   0/1/2 — index into the IR-channel axis.
    tict : np.ndarray  Smoothed ICT temperature from :func:`extract_ir_telemetry`.
    mask : np.ndarray  Boolean scanline mask.
    ict_threshold, space_threshold : float
    window : int

    Returns
    -------
    cs : np.ndarray  Smoothed space counts (scanlines,).
    cict : np.ndarray  Smoothed ICT counts (scanlines,).
    ce : np.ndarray  Earth counts (scanlines, pixels).
    """
    space = ds["full_space_counts"].isel(channel_name=(cal_index - 3)).mean(axis=1).values
    ict = ds["full_ict_counts"].isel(ir_channel_name=cal_index).mean(axis=1).values
    ce = ds["counts"].values[:, :, cal_index - 3]

    gd = ~np.isfinite(space)
    if np.sum(gd) > 0:
        space[gd] = 0
    gd = ~np.isfinite(ict)
    if np.sum(gd) > 0:
        ict[gd] = 0

    ict[mask] = 0
    space[mask] = 0

    zeros = ict < ict_threshold
    nonzeros = ~zeros
    no37 = False
    try:
        ict[zeros] = np.interp(zeros.nonzero()[0], nonzeros.nonzero()[0], ict[nonzeros])
    except ValueError:
        no37 = True

    if not no37:
        zeros = space < space_threshold
        nonzeros = ~zeros
        space[zeros] = np.interp(zeros.nonzero()[0], nonzeros.nonzero()[0], space[nonzeros])
    else:
        space[:] = np.nan
        ict[:] = np.nan

    weighting_function = np.ones(window, dtype=float) / window
    half = (window - 1) // 2
    half_up = (window + 1) // 2

    def _convolve_clamp(arr):
        out = np.convolve(arr, weighting_function, "same")
        out[:half] = out[half]
        out[-half:] = out[-half_up]
        return out

    return _convolve_clamp(space), _convolve_clamp(ict), ce


def allan_deviation(space, bad_scan=None):
    """Determine the Allan deviation (noise) from space view counts filtering
    out bad space view lines. Written by J.Mittaz, University of Reading"""

    if len(space.shape) != 2:
        raise Exception("utils.allan_deviation input space view not 2-dimensional")
    #
    # Get good scanlines if filter present
    #
    if bad_scan is not None:
        gd = (bad_scan == 0)
        newsp = space[gd,:]
    else:
        newsp = space

    #
    # Allan deviation is sqrt of allan variance which is
    #
    #        allan_variance = 0.5 <(y_n+1-y_n)**2>
    #
    allan_variance = 0.5 * np.mean(np.diff(newsp, axis=1) ** 2)
    #
    # Return the Allan deviation in counts (sqrt Allan variance)
    #
    return np.sqrt(allan_variance)

class convBT:
    """Routine to convert temperature to radiance and vice-versa."""
    def t_to_rad(self,tprt):

        tsBB = self.A + self.B*tprt
        return self.nBB_num / (np.exp(self.c2_nu_c / tsBB) - 1.0)

    def rad_to_t(self,rad):

        corrT = self.c2_nu_c/np.log((self.nBB_num/rad)+1.)
        return (corrT-self.A)/self.B

    def rad_to_t_uncert(self,rad,urad):

        T = self.rad_to_t(rad)
        T1 = self.rad_to_t(rad+urad)
        T2 = self.rad_to_t(rad-urad)

        return T,(np.abs(T-T1)+np.abs(T-T2))/2.

    def t_to_rad_uncert(self,T,uT):

        rad = self.t_to_rad(T)
        rad1 = self.t_to_rad(T+uT)
        rad2 = self.t_to_rad(T-uT)

        return rad,(np.abs(rad-rad1)+np.abs(rad-rad2))/2.

    def __init__(self,cal,chan):

        # constants
        self.c1 = 1.1910427e-5  # mW/m^2/sr/cm^{-4}
        self.c2 = 1.4387752  # cm K
        # coefficients
        self.A = cal.to_eff_blackbody_intercept[chan]
        self.B = cal.to_eff_blackbody_slope[chan]
        self.nu_c = cal.centroid_wavenumber[chan]
        self.nBB_num = self.c1 * (self.nu_c**3)
        self.c2_nu_c = self.c2 * self.nu_c

def get_bad_space_counts(sp_data, ict_data=None):
    """Find bad space count data.

    Space count data is voltage clamped so should have very close to the same value close to 950 - 960.
    """

    #
    # Use robust estimators to get thresholds for space counts
    # Use 4 sigma threshold from median
    # Ensure only for good data
    #
    gd = np.isfinite(sp_data)
    if np.sum(gd) == 0:
        sp_bad_data = np.zeros(sp_data.shape,dtype=bool)
        sp_bad_data[:,:] = True
        return sp_bad_data

    quantile = np.quantile(sp_data[gd].flatten(),[0.25,0.75])
    if quantile[0] == quantile[1]:
        quantile[1] = quantile[1]+0.5
    std = (quantile[1]-quantile[0])/1.349
    sp_bad_data = np.zeros(sp_data.shape,dtype=bool)
    sp_bad_data[:,:] = True
    gd = np.isfinite(sp_data)
    if ict_data is not None:
        sp_bad_data[gd] = ~((np.abs(sp_data[gd] - np.median(sp_data[gd].flatten()))/\
                             std < 5.)&(ict_data[gd] > 0))
    else:
        sp_bad_data[gd] = ~((np.abs(sp_data[gd] - np.median(sp_data[gd].flatten()))/\
                             std < 5.))

    return sp_bad_data


def get_noise(total_space, total_ict, window, twelve_micron):
    """Get noise estimates from the counts.

    .. deprecated::
        Use :func:`_compute_noise_arrays` with a list of :class:`IRChannelSpec`
        objects instead.  This wrapper exists for backward compatibility only
        and will be removed in a future release.
    """
    import warnings
    from types import SimpleNamespace
    warnings.warn(
        "get_noise() is deprecated. Use _compute_noise_arrays() instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    specs = [SimpleNamespace(cal_index=0), SimpleNamespace(cal_index=1)]
    if twelve_micron:
        specs.append(SimpleNamespace(cal_index=2))

    bad_scans = _compute_bad_scans(specs, total_space, total_ict)

    def _noise_triple(cal_idx):
        return _channel_noise(
            total_space[:, :, cal_idx],
            total_ict[:, :, cal_idx],
            bad_scans,
            window,
        )

    n1, av1, av_ict1 = _noise_triple(0)
    n2, av2, av_ict2 = _noise_triple(1)
    n3 = av3 = av_ict3 = None
    if twelve_micron:
        n3, av3, av_ict3 = _noise_triple(2)

    return n1, n2, n3, av1, av2, av3, av_ict1, av_ict2, av_ict3, bad_scans

def smooth_data(y,length):
    """Smooth data over given length"""

    outy = np.zeros(len(y),dtype=y.dtype)
    leny = len(y)-1
    for i in range(leny+1):
        minx=max([i-length,0])
        maxx=min([i+length,leny])
        outy[i] = np.mean(y[minx:maxx+1])

    return outy

def get_uICT(gainval,CS,CICT,Tict,NS,convT,bad_scans,solar_scans,window):
    """Get ICT temperature gradient uncertainty based on analysis of the
    gain variations in the 3.7 micron channel
    Estimate on "window" scanline length
    Only uses 'good' data"""

    #
    # Only use good data
    #
    gd = (bad_scans == 0)&(solar_scans == 0)
    Lict = convT.t_to_rad(Tict[gd])

    gain = (Lict-NS)/(CS[gd]-CICT[gd])
    sp_ict = CS[gd]-CICT[gd]
    Tcorr = convT.A+convT.B*Tict[gd]

    dGain_dICT = (convT.B/sp_ict)*\
        convT.nBB_num*np.exp(convT.c2_nu_c/Tcorr)*\
        (convT.c2_nu_c/(Tcorr**2))/\
        (np.exp(convT.c2_nu_c/Tcorr)-1.)**2

    dgain = (gain-gainval)

    dT = np.zeros(len(CS))
    dT[:] = np.nan
    #
    # Calculate delta ICT for 'good' cases
    #
    dT[gd] = dgain/dGain_dICT

    uICT = np.zeros(len(CS))
    for i in range(len(dT)):
        if np.isfinite(dT[i]):
            minx = max([0,i-window//2])
            maxx = min([len(dT)+1,i+window//2+1])
            deltaT = dT[minx:maxx]
            gd = np.isfinite(deltaT)
            if np.sum(gd) > 0:
                uICT[i] = np.std(deltaT[gd])
            else:
                uICT[i] = np.nan
        else:
            uICT[i] = np.nan

    return uICT

def get_ict_uncert(Tict, prt_random, prt_bias, uICT, convT):
    """Get uncertainty in radiance/temperature of ICT on the basis of ICT uncertainties.

    Accepts scalar or array inputs for *Tict* and *uICT*.
    """
    # Note: in operational calibration PRTs are averaged over 4,
    # so input PRT uncertainties need to be divided by 2.
    rad, urand = convT.t_to_rad_uncert(Tict, prt_random / 2.)

    # Systematic does not average down; NaN uICT → NaN prt_sys.
    # Works element-wise for both scalars and arrays.
    fin = np.isfinite(uICT)
    safe_u = np.where(fin, uICT, 0.)
    prt_sys = np.where(fin, np.sqrt(prt_bias**2 + safe_u**2), np.nan)

    return urand, prt_sys, rad

class fit_ict_pars(object):
    """Model for ICT gradients"""

    def model(self,X,a0,a1,a2,a3):

        ICT = X[:,0] + a0*X[:,1] + a1*X[:,2] + a2*X[:,3] + a3*X[:,4]
        rad_ict = self.convT.t_to_rad(ICT)

        return (rad_ict - X[:,7])/(X[:,5] - X[:,6])

    def __init__(self,convT):

        self.convT = convT

def find_solar_ind(position,new_gain2,new_cs,new_ct,new_ict,new_ict1,
                   new_ict2,new_ict3,new_ict4,solZA,first=True):
    """Find possible solar contamination on one side of the solar max
    position"""

    pos = np.nonzero(solZA == solZA.max())[0][0]
    #
    # Look before or after max solZA location
    #
    if first:
        posmin = np.nonzero(solZA[0:pos] == solZA[0:pos].min())[0][0]
    else:
        posmin = pos+np.nonzero(solZA[pos:] == solZA[pos:].min())[0][0]
    if pos < posmin:
        posn = position[pos:posmin]
        newg = new_gain2[pos:posmin]
        sza = solZA[pos:posmin]
    else:
        posn = position[posmin:pos]
        newg = new_gain2[posmin:pos]
        sza = solZA[posmin:pos]

    #
    # Only look in range satellite ZA >= 100 and <= 125
    #
    gd = (sza >= 100)&(sza <= 125)&np.isfinite(newg)
    if np.sum(gd) == 0:
        print("ERROR: No data available to find solar contamination")
        print("       MinSZA = {0:f} MaxSZA = {1:f}".format(sza.min(),sza.max()))
        return None,None,None
    posn2 = posn[gd]
    pos2 = np.nonzero(newg[gd] == newg[gd].min())[0][0]
    newg = newg[gd]
    peak_gain = newg[pos2]
    #
    # Track back to a defined limit when gradient is reversed or zero
    #
    peak_location = posn2[pos2]
    ok1 = False

    window_size = 40
    while not ok1:
        side1 = 0
        for i in range(pos2,0,-1):
            if i > pos2-window_size:
                continue
            minx = max([0,i-window_size])
            maxx = min([peak_location,i+window_size])
            y = newg[minx:maxx+1]
            x = np.arange(len(y))
            p = np.polyfit(x,y,1)
            if p[0] >= 0.:
                side1 = i
                break
        #
        # Check that we are not at an edge
        #
        if side1 == 0:
            print("Cannot find peak for solar contamination")
            return None,None,None

        if pos2+window_size < len(newg):
            ok1 = True
        else:
            #
            # Find next peak < side1 position
            #
            pos2 = np.nonzero(newg[:side1] == newg[:side1].min())[0][0]
            peak_location = posn2[pos2]

    #
    # Look at other side - shouldn't have same problem as side1 with possible
    # erroneous peak
    #
    side2 = len(newg)-1
    for i in range(pos2,len(newg)):
        if i < pos2+window_size:
            continue
        minx = max([0,i-window_size])
        maxx = min([peak_location,i+window_size])
        y = newg[minx:maxx+1]
        x = np.arange(len(y))
        p = np.polyfit(x,y,1)
        if p[0] <= 0.:
            side2 = i
            break

    #
    # Add Take larger difference case
    #
    diff1 = np.abs(side1-pos2)
    diff2 = np.abs(side2-pos2)
    if diff1 > diff2:
        side1 = diff1
        side2 = diff1
    else:
        side1 = diff2
        side2 = diff2

    #
    # Check to see if we are already at a background value or have to go
    # further (case where second peak on side of signal
    #
    newpos1 = max([0,pos2-side1])
    newpos2 = min([len(newg)-1,pos2+side2])
    maxdiff1 = np.abs(peak_gain-newg[newpos1])
    maxdiff2 = np.abs(peak_gain-newg[newpos2])
    maxdiff = max([maxdiff1,maxdiff2])
    update_side1 = False
    for i in range(side1+side1//2):
        bloc = pos2-side1-i
        if bloc <= 0:
            break
        test_ratio = maxdiff/np.abs(peak_gain-newg[bloc])
        if test_ratio > 0.7 and test_ratio < 1.3:
            update_side1 = True
            side1 = bloc
            break

    update_side2 = False
    for i in range(side2+side2//2):
        bloc = pos2+side2+i
        if bloc >= len(newg):
            break
        test_ratio = maxdiff/np.abs(peak_gain-newg[bloc])
        if test_ratio > 0.7 and test_ratio < 1.3:
            update_side2 = True
            side2 = bloc
            break

    #
    # Reset to full location
    #
    if not update_side1:
        side1 = newpos1
    if not update_side2:
        side2 = newpos2

    #
    # Can't get background so eject
    #
    if side1 == 0:
        print("ERROR: Cannot get a background for solar contamination")
        return None,None,None

    #
    # Background and sigma
    #
    maxx = side1
    minx = max([0,maxx-300])
    backdata1 = newg[minx:maxx]
    pos_1 = (minx+maxx)/2.
    gd = np.isfinite(backdata1)
    if np.sum(gd) == 0:
        print("ERROR: Cannot get background for solar contamination")
        return None,None,None
    backg1 = np.mean(backdata1[gd])
    backg1_std = np.std(backdata1[gd])

    maxx = side2
    minx = min([len(newg)-1,maxx+300])
    backdata2 = newg[maxx:minx]
    gd = np.isfinite(backdata2)
    if np.sum(gd) == 0:
        pos_2 = -1
    else:
        backg2 = np.mean(backdata2[gd])
        pos_2 = (minx+maxx)/2.
        backg2_std = np.std(backdata2[gd])

    #
    # If pos_2 there, and OK interpolate
    #
    if pos_2 > -1:
        #
        # Make sure back2 is < 50% of peak
        #
        if (backg1-peak_gain)*0.25 > (backg2-peak_gain):
            #
            # Less than 25% from peak so dont continue
            #
            print("ERROR: back2 is too high (solar contamination)")
            return None,None,None
        slope = (backg2-backg1)/(pos_2-pos_1)
        cnst = backg1 - slope*pos_1
        backg = cnst + slope*pos2
        backg_std = (backg1_std+backg2_std)/2.
    else:
        print("ERROR: No back2 value (solar contamination)")
        return None,None,None

    sigma = np.abs(peak_gain-backg)/backg_std
    if sigma <= 7.5:
        return None,None,None

    side1 = posn2[side1]
    side2 = posn2[side2]

    return side1,side2,peak_location

def find_solar(ds,mask,convT=None):
    """Find solar contamination/variable gain points for GAC and for ~full
    orbits"""

    #
    # Get parameters for kernals/uncertainties
    #
    window, prt_bias, prt_sys, prt_threshold, ict_threshold, \
        space_threshold = get_uncert_parameter_thresholds()
    if ds["channels"].values.shape[1] == 409:
        gacdata = True
    else:
        gacdata = False
        #
        # Code only works for GAC
        #
        print("ERROR: Solar contamination detection only works for GAC data")
        return None,None,None,None

    #
    # Get calibration coefficients for 3.7 micron channel
    #
    cal = Calibrator(
        ds.attrs["spacecraft_name"])
    NS_1 = cal.space_radiance[0]

    if convT is None:
        convT = convBT(cal,0)

    CS_1,CICT_1,CE_1,Tict,ict1,ict2,ict3,ict4,solZAin \
        = get_vars(ds,0,convT,
                   window,
                   prt_threshold,
                   ict_threshold,
                   space_threshold,
                   gacdata,
                   cal,
                   mask,
                   out_prt=True,
                   out_solza=True)

    meanT = (ict1+ict2+ict3+ict4)/4.
    radBB = convT.t_to_rad(meanT)
    gain3 = (radBB-NS_1)/(CS_1-CICT_1)

    #
    # Get location of min PRT stdev as closest to zero error case
    #
    X = np.zeros((len(meanT),4))
    X[:,0] = ict1
    X[:,1] = ict2
    X[:,2] = ict3
    X[:,3] = ict4
    stdev = np.std(X,axis=1)
    try:
        minstd_pos = np.nonzero(stdev == stdev.min())[0][0]
    except IndexError:
        minstd_pos = np.nonzero(stdev == stdev.min())[0]

    #
    # Get PRT differences from Mean
    #
    prt_diff1 = ict1 - meanT
    prt_diff2 = ict2 - meanT
    prt_diff3 = ict3 - meanT
    prt_diff4 = ict4 - meanT

    #
    # Now fit model
    #
    X = np.zeros((len(meanT),8))
    Y = np.zeros((len(meanT)))
    X[:,0] = meanT[:]
    X[:,1] = prt_diff1[:]
    X[:,2] = prt_diff2[:]
    X[:,3] = prt_diff3[:]
    X[:,4] = prt_diff4[:]
    X[:,5] = CS_1[:]
    X[:,6] = CICT_1[:]
    X[:,7] = NS_1
    Y[:] = gain3[minstd_pos]

    model = fit_ict_pars(convT)
    p,covar = curve_fit(model.model,X,Y,p0=[0.,0.,0.,0.])

    #
    # Get updated gain
    #
    new_gain = model.model(X,p[0],p[1],p[2],p[3])

    #
    # Find peak gain within solZA 100-125 degrees
    #
    # First find range between min/max values
    #
    position = np.arange(len(solZAin)).astype(dtype=np.int32)
    gd = np.isfinite(solZAin)
    position = position[gd]
    new_gain2 = new_gain[gd]
    new_cs = CS_1[gd]
    new_ct = CICT_1[gd]
    new_ict = Tict[gd]
    new_ict1 = ict1[gd]
    new_ict2 = ict2[gd]
    new_ict3 = ict3[gd]
    new_ict4 = ict4[gd]
    solZA = solZAin[gd]

    #
    # Now look before and after solZA max position
    #
    #
    # One side of solZA max
    #
    side1_1,side2_1,peak_location_1 = find_solar_ind(position,new_gain2,
                                                     new_cs,new_ct,new_ict,
                                                     new_ict1,new_ict2,
                                                     new_ict3,new_ict4,
                                                     solZA,first=True)
    #
    # Other side of solZA max
    #
    side1_2,side2_2,peak_location_2 = find_solar_ind(position,new_gain2,
                                                     new_cs,new_ct,new_ict,
                                                     new_ict1,new_ict2,
                                                     new_ict3,new_ict4,
                                                     solZA,first=False)

    if side1_1 is not None and side1_2 is not None:
        return side1_1,side2_1,peak_location_1,\
            solZAin[peak_location_1],side1_2,side2_2,\
            peak_location_2,solZAin[peak_location_2]
    elif side1_1 is not None and side1_2 is None:
        return side1_1,side2_1,peak_location_1,\
            solZAin[peak_location_1],None,None,None,None
    elif side1_1 is None and side1_2 is not None:
        return None,None,None,None,side1_2,side2_2,\
            peak_location_2,solZAin[peak_location_2]
    else:
        return None,None,None,None,None,None,None,None


def _radiance_random_uncert(spec, noise, av_noise, ict_noise, ict_random,
                             Lict, CS, CE, CICT):
    """Random radiance uncertainty for one IR channel.

    Replaces :func:`get_random` — uses :class:`IRChannelSpec` instead of a
    magic ``channel`` integer, so the 3.7µm / 11-12µm branching is expressed
    as data (``spec.has_nonlinear``) rather than control flow.

    Works on scalars or NumPy arrays; broadcasting follows normal NumPy rules.
    CS, CICT must be broadcastable with CE.
    """
    NS = spec.space_radiance
    _, c1, c2 = spec.nonlin_coeffs

    dLlin_dCS = (Lict - NS) * (CS - CE) / (CS - CICT) ** 2 + (Lict - NS) / (CS - CICT)
    dLlin_dCICT = -(Lict - NS) * (CS - CE) / (CS - CICT) ** 2
    dLlin_dCE = -(Lict - NS) / (CS - CICT)
    dLlin_dLict = (CS - CE) / (CS - CICT)

    if spec.has_nonlinear:
        Llin = NS + (Lict - NS) * (CS - CE) / (CS - CICT)
        scale = 1.0 + c1 + c2 * Llin
        dL_dCS = dLlin_dCS * scale
        dL_dCICT = dLlin_dCICT * scale
        dL_dCE = dLlin_dCE * scale
        dL_dLict = dLlin_dLict * scale
    else:
        dL_dCS, dL_dCICT, dL_dCE, dL_dLict = dLlin_dCS, dLlin_dCICT, dLlin_dCE, dLlin_dLict

    uncert = (
        (dL_dCS ** 2) * (av_noise ** 2)
        + (dL_dCICT ** 2) * (ict_noise ** 2)
        + (dL_dCE ** 2) * (noise ** 2)
        + (dL_dLict ** 2) * (ict_random ** 2)
    )
    return np.sqrt(uncert)


def _radiance_sys_uncert(spec, uICT, Tict, CS, CE, CICT):
    """Systematic radiance uncertainty for one IR channel.

    Replaces :func:`get_sys` — uses :class:`IRChannelSpec` instead of a magic
    ``channel`` integer.  Returns NaN when *uICT* is NaN (no valid systematic
    uncertainty); the ``sys_there`` boolean flag from the old API is gone —
    callers detect the NaN case with ``np.isfinite``.
    """
    NS = spec.space_radiance
    _, c1, c2 = spec.nonlin_coeffs

    if not np.isfinite(uICT).any():
        return np.full(CE.shape, np.nan)

    Lict, uradTict = spec.conv.t_to_rad_uncert(Tict, uICT)
    dLlin_dLict = (CS - CE) / (CS - CICT)

    if spec.has_nonlinear:
        Llin = NS + (Lict - NS) * (CS - CE) / (CS - CICT)
        dL_dLict = dLlin_dLict * (1.0 + c1 + c2 * Llin)
    else:
        dL_dLict = dLlin_dLict

    return np.sqrt((dL_dLict ** 2) * (uradTict ** 2))


def get_random(channel, noise, av_noise, ict_noise, ict_random, Lict, CS, CE, CICT, NS,
               c1, c2):
    """Deprecated — use :func:`_radiance_random_uncert` instead."""
    import warnings
    warnings.warn(
        "get_random is deprecated; use _radiance_random_uncert(spec, ...) instead.",
        DeprecationWarning, stacklevel=2,
    )
    from types import SimpleNamespace
    spec = SimpleNamespace(
        space_radiance=NS,
        nonlin_coeffs=(0., c1, c2),
        has_nonlinear=(channel != 1),
    )
    return _radiance_random_uncert(spec, noise, av_noise, ict_noise, ict_random,
                                   Lict, CS, CE, CICT)

def get_sys(channel, uICT, Tict, CS, CE, CICT, NS, c1, c2, convT):
    """Deprecated — use :func:`_radiance_sys_uncert` instead."""
    import warnings
    warnings.warn(
        "get_sys is deprecated; use _radiance_sys_uncert(spec, ...) instead.",
        DeprecationWarning, stacklevel=2,
    )
    from types import SimpleNamespace
    spec = SimpleNamespace(
        space_radiance=NS,
        nonlin_coeffs=(0., c1, c2),
        has_nonlinear=(channel != 1),
        conv=convT,
    )
    result = _radiance_sys_uncert(spec, uICT, Tict, CS, CE, CICT)
    sys_there = np.isfinite(result).any()
    return result, sys_there

def get_vars(ds,channel,convT,wlength,prt_threshold,ict_threshold,
             space_threshold,gac,cal,mask,out_prt=False,out_solza=False):
    """Get variables from xarray including smoothing and interpolation"""
    space = ds["full_space_counts"].isel(channel_name=(channel - 3)).mean(axis=1).values
    # Defensive copy: .values[:] returns a *view* and we mutate prt below
    # (prt[gd] = 0, prt[ifix] = np.interp(...)). Without the copy these
    # writes propagate back into ds["mean_prt_counts"], which (a) corrupts
    # the dataset for any later caller and (b) makes the 2nd/3rd
    # per-channel get_vars() call see different inputs from the 1st.
    prt = ds["mean_prt_counts"].values.copy()
    ict = ds["full_ict_counts"].isel(ir_channel_name=channel).mean(axis=1).values
    ce = ds["counts"].values[:,:,channel - 3]
    midpoint = ds["sun_zen"].shape[1]//2
    line_numbers = ds["scan_line_index"].data

    if out_solza:
        solza = ds["sun_zen"].values[:,midpoint]

    #
    # Set nan's to value to be caught by interpolation routines
    #
    gd = ~np.isfinite(prt)
    if np.sum(gd) > 0:
        prt[gd] = 0
    gd = ~np.isfinite(ict)
    if np.sum(gd) > 0:
        ict[gd] = 0
    gd = ~np.isfinite(space)
    if np.sum(gd) > 0:
        space[gd] = 0

    #
    # Removed old way of doing PRT indexing
    #
    # PRT index check
    #
    # PRTs. See reader.get_telemetry implementations.
    #
    #for offset in range(5):
    #    # According to the KLM Guide the fill value between PRT measurments is 0, but we search
    #    # for the first measurement gap using the threshold, because the fill value is in practice
    #    # not always exactly 0.
    #    if np.median(prt[(line_numbers - line_numbers[0]) % 5 == offset]) < prt_threshold:
    #        break
    #    else:
    #        raise IndexError("No PRT 0-index found!")
    #
    # get the PRT index, iprt equals to 0 corresponds to the measurement gaps
    # This can give wrong ICT temperatures
    #
    #iprt = (line_numbers - line_numbers[0] + 5 - offset) % 5

    # Get PRT mapping using new technique
    iprt = get_prt_nos(prt,prt_threshold,line_numbers,gac)

    #
    # Interpolate over bad prt values - from pygac calibrate_thermal
    #
    # fill measured values below threshold by interpolation
    #
    ifix = np.where(np.logical_and(iprt == 1, prt <= prt_threshold))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 1, prt > prt_threshold))
        if len(inofix[0]):
            prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])
        else:
            raise IndexError("No good prt1 data")

    ifix = np.where(np.logical_and(iprt == 2, prt <= prt_threshold))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 2, prt > prt_threshold))
        if len(inofix[0]):
            prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])
        else:
            raise IndexError("No good prt2 data")

    ifix = np.where(np.logical_and(iprt == 3, prt <= prt_threshold))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 3, prt > prt_threshold))
        if len(inofix[0]):
            prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])
        else:
            raise IndexError("No good prt3 data")

    ifix = np.where(np.logical_and(iprt == 4, prt <= prt_threshold))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 4, prt > prt_threshold))
        if len(inofix[0]):
            prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])
        else:
            raise IndexError("No good prt4 data")

    #
    # Convert to temperature
    #
    # calculate PRT temperature using equation (7.1.2.4-1) KLM Guide
    # Tprt = d0 + d1*Cprt + d2*Cprt^2 + d3*Cprt^3 + d4*Cprt^4
    # Note: First dimension of cal.d are the five coefficient indicees
    #
    tprt = np.polynomial.polynomial.polyval(prt, cal.d[:, iprt], tensor=False)

    #
    # Get interpolated values as done in pygac
    #
    tprt_interp = np.copy(tprt)
    zeros = iprt == 0
    nonzeros = np.logical_not(zeros)

    tprt_interp[zeros] = np.interp((zeros).nonzero()[0],
                            (nonzeros).nonzero()[0],
                            tprt[nonzeros])
    #
    # Interpolate over each PRT number
    #
    tprt1_interp = np.copy(tprt)
    zeros = (iprt == 0)|(iprt != 1)
    nonzeros = np.logical_not(zeros)

    tprt1_interp[zeros] = np.interp((zeros).nonzero()[0],
                            (nonzeros).nonzero()[0],
                            tprt[nonzeros])

    tprt2_interp = np.copy(tprt)
    zeros = (iprt == 0)|(iprt != 2)
    nonzeros = np.logical_not(zeros)

    tprt2_interp[zeros] = np.interp((zeros).nonzero()[0],
                            (nonzeros).nonzero()[0],
                            tprt[nonzeros])

    tprt3_interp = np.copy(tprt)
    zeros = (iprt == 0)|(iprt != 3)
    nonzeros = np.logical_not(zeros)

    tprt3_interp[zeros] = np.interp((zeros).nonzero()[0],
                            (nonzeros).nonzero()[0],
                            tprt[nonzeros])

    tprt4_interp = np.copy(tprt)
    zeros = (iprt == 0)|(iprt != 4)
    nonzeros = np.logical_not(zeros)

    tprt4_interp[zeros] = np.interp((zeros).nonzero()[0],
                            (nonzeros).nonzero()[0],
                            tprt[nonzeros])

    # Thresholds to flag missing/wrong data for interpolation
    # Remove masked data
    ict[mask] = 0
    space[mask] = 0
    zeros = ict < ict_threshold
    nonzeros = np.logical_not(zeros)
    no37 = False
    try:
        ict[zeros] = np.interp((zeros).nonzero()[0],
                               (nonzeros).nonzero()[0],
                               ict[nonzeros])
    except ValueError: # 3b has no valid data
        no37 = True
    if not no37:
        zeros = space < space_threshold
        nonzeros = np.logical_not(zeros)

        space[zeros] = np.interp((zeros).nonzero()[0],
                                 (nonzeros).nonzero()[0],
                                 space[nonzeros])
    else:
        space[:] = np.nan
        ict[:] = np.nan

    #
    # Make averages and do using pygacs method at this point
    #
    weighting_function = np.ones(wlength, dtype=float) / wlength
    tprt_convolved = np.convolve(tprt_interp, weighting_function, "same")
    tprt1_convolved = np.convolve(tprt1_interp, weighting_function, "same")
    tprt2_convolved = np.convolve(tprt2_interp, weighting_function, "same")
    tprt3_convolved = np.convolve(tprt3_interp, weighting_function, "same")
    tprt4_convolved = np.convolve(tprt4_interp, weighting_function, "same")
    ict_convolved = np.convolve(ict, weighting_function, "same")
    space_convolved = np.convolve(space, weighting_function, "same")

    # take care of the beginning and end
    tprt_convolved[0:(wlength - 1) // 2] = tprt_convolved[(wlength - 1) // 2]
    tprt1_convolved[0:(wlength - 1) // 2] = tprt1_convolved[(wlength - 1) // 2]
    tprt2_convolved[0:(wlength - 1) // 2] = tprt2_convolved[(wlength - 1) // 2]
    tprt3_convolved[0:(wlength - 1) // 2] = tprt3_convolved[(wlength - 1) // 2]
    tprt4_convolved[0:(wlength - 1) // 2] = tprt4_convolved[(wlength - 1) // 2]
    ict_convolved[0:(wlength - 1) // 2] = ict_convolved[(wlength - 1) // 2]
    space_convolved[0:(wlength - 1) // 2] = space_convolved[(wlength - 1) // 2]
    tprt_convolved[-(wlength - 1) // 2:] = tprt_convolved[-((wlength + 1) // 2)]
    tprt1_convolved[-(wlength - 1) // 2:] = tprt1_convolved[-((wlength + 1) // 2)]
    tprt2_convolved[-(wlength - 1) // 2:] = tprt2_convolved[-((wlength + 1) // 2)]
    tprt3_convolved[-(wlength - 1) // 2:] = tprt3_convolved[-((wlength + 1) // 2)]
    tprt4_convolved[-(wlength - 1) // 2:] = tprt4_convolved[-((wlength + 1) // 2)]
    ict_convolved[-(wlength - 1) // 2:] = ict_convolved[-((wlength + 1) // 2)]
    space_convolved[-(wlength - 1) // 2:] = space_convolved[-((wlength + 1) // 2)]

    if out_prt:
        if out_solza:
            return space_convolved,ict_convolved,ce,tprt_convolved,\
                tprt1_convolved,tprt2_convolved,tprt3_convolved,\
                tprt4_convolved,solza
        else:
            return space_convolved,ict_convolved,ce,tprt_convolved,\
                tprt1_convolved,tprt2_convolved,tprt3_convolved,\
                tprt4_convolved
    else:
        if out_solza:
            return space_convolved,ict_convolved,ce,tprt_convolved,solza
        else:
            return space_convolved,ict_convolved,ce,tprt_convolved

@contextmanager
def open_zenodo_uncert_file(platform, decode_times=True):
    import ssl

    import fsspec
    import truststore

    #
    # Force the right naming conventions
    #
    if platform == "noaa6":
        platform = "noaa06"
    elif platform == "noaa7":
        platform = "noaa07"
    elif platform == "noaa8":
        platform = "noaa08"
    elif platform == "noaa9":
        platform = "noaa09"

    ctx = truststore.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
    coef_file = fsspec.open_local(f"simplecache::https://zenodo.org/records/16926055/files/{platform}_uncert.nc#mode=bytes",
                                  simplecache=dict(cache_storage=gettempdir(), same_names=True),
                                  https=dict(ssl=ctx))
    with xr.open_dataset(coef_file, decode_times=decode_times) as d:
        yield d


def get_gainval(time,intimes,avhrr,prt1,prt2,prt3,prt4,CS,CICT,CE,NS,
                bad_scan,convT,window,calculate=False):
    """Estimate gain value at smallest stdev point in orbit either from
    file or estimate it from data"""
    #
    # If no 3.7 micron data present then don't do anything
    #
    if np.sum(np.isfinite(CS)) == 0:
        print("ERROR: No 3.7 micron data for gain calc. so no ICT uncertainty")
        return None,None

    #
    # Use nearest in time if HRPT (not calculate)
    #
    if not calculate:
        #
        # Open file containing 3.7mu gain value and interpolate over time
        #
        try:
            with open_zenodo_uncert_file(avhrr) as d:
                intime = d["time_gain"].values[:]
                ingain = d["gain"].values[:]
        except FileNotFoundError:
            raise Exception("ERROR: Gain can not be determined because zenodo not available")

        if time > intime[-1]:
            print("ERROR: file time > last time with max gain values (HRPR/LAC)")
            return None,None

        timediff = (intime-time)/np.timedelta64(1,"s")
        timediff = np.abs(timediff)
        timediff_min = timediff.min()
        #
        # HRPT/LAC
        #
        pos = np.nonzero(timediff == timediff_min)[0][0]
        return ingain[pos],intime[pos]
    else:
        #
        # No nearby gain estimate or force calculate so calculate
        # from data using min std of prts
        #
        gd = (bad_scan == 0)
        times = intimes[gd]
        prt1 = prt1[gd]
        prt2 = prt2[gd]
        prt3 = prt3[gd]
        prt4 = prt4[gd]
        CS = CS[gd]
        CICT = CICT[gd]
        #
        # Get stdev
        #
        X = np.zeros((len(prt1),4))
        X[:,0] = prt1
        X[:,1] = prt2
        X[:,2] = prt3
        X[:,3] = prt4
        stdev = np.std(X,axis=1)
        pos = np.nonzero(stdev == stdev.min())[0][0]
        Lict = convT.t_to_rad(np.mean(X,axis=1))
        gain = (Lict-NS)/(CS-CICT)
        return gain[pos],times[pos]

def get_pixel(Lict,CS,CE,CICT,NS,c0,c1,c2):
    """Get radiance of pixel using calibration"""

    Llin = NS + (Lict-NS)*(CS-CE)/(CS-CICT)
    if NS != 0.:
        LE = Llin + c0 + c1*Llin + c2*Llin*Llin
        return LE
    else:
        return Llin

def get_uncert_parameter_thresholds(vischans=False):
    """Return required constants for IR and visible uncertainty cases. Single
    point so any changes will be correctly applied across both"""
    #
    # Define averaging kernel based on value in noaa.py
    # Also set PRT uncertainty components and thresholds
    #
    if vischans:
        window = 51
        solar_contam_threshold = 0.05
        solar_contam_sza_threshold = 102.
        return window,solar_contam_threshold,solar_contam_sza_threshold
    else:
        window = 51
        prt_bias = 0.01
        prt_sys = 0.1
        prt_threshold = 50
        ict_threshold = 100
        space_threshold = 100

        return window, prt_bias, prt_sys, prt_threshold, ict_threshold, \
            space_threshold

def get_solar_from_file(platform, ds):
    """Read in possible solar contamination times from uncertainty files.

    Only used for LAC/HRPT data where there is not enough data to detect
    possible solar contamination so GAC estimates are used.
    """

    #
    # Get times in seconds from
    #
    time = (ds["times"].values[:] - \
            np.datetime64("1970-01-01T00:00:00"))/\
            np.timedelta64(1,"s")

    #
    # Read file
    #
    try:
        with open_zenodo_uncert_file(platform, decode_times=False) as d:
            solar_start_time_1 = d["gain1_solar_start"].values[:]
            solar_stop_time_1 = d["gain1_solar_stop"].values[:]
            solar_start_time_2 = d["gain2_solar_start"].values[:]
            solar_stop_time_2 = d["gain2_solar_stop"].values[:]
    except FileNotFoundError:
        raise Exception("ERROR: Solar data can mot be determined because zenodo not available")

    gd = np.isfinite(solar_start_time_1)&np.isfinite(solar_stop_time_1)
    solar_start_time_1 = solar_start_time_1[gd]
    solar_stop_time_1 = solar_stop_time_1[gd]
    gd = np.isfinite(solar_start_time_2)&np.isfinite(solar_stop_time_2)
    solar_start_time_2 = solar_start_time_2[gd]
    solar_stop_time_2 = solar_stop_time_2[gd]
    #
    # Match to times in file
    #
    min_solar_1 = -1
    max_solar_1 = -1
    for i in range(len(solar_start_time_1)):
        gd = (time >= solar_start_time_1[i])&(time <= solar_stop_time_1[i])
        if np.sum(gd) > 0:
            index = np.arange(len(time)).astype(dtype=np.int32)
            index = index[gd]
            min_solar_1 = index[0]
            max_solar_1 = index[-1]
    min_solar_2 = -1
    max_solar_2 = -1
    for i in range(len(solar_start_time_2)):
        gd = (time >= solar_start_time_2[i])&(time <= solar_stop_time_2[i])
        if np.sum(gd) > 0:
            index = np.arange(len(time)).astype(dtype=np.int32)
            index = index[gd]
            min_solar_2 = index[0]
            max_solar_2 = index[-1]

    #
    # Return solar location if available
    #
    return min_solar_1, max_solar_1, min_solar_2, max_solar_2

def ir_uncertainty(ds,mask):
    """Create the uncertainty components for the IR channels. These include

    1) Random
          a) Noise
          b) Digitisation
          c) ICT PRT Noise
    2) Systematic
          a) ICT Temperature uncertainty
          b) PRT Bias
          c) Calibration coefs/measurement equation uncertainty

    Inputs:
            ds : Input xarray dataset containing data for calibration
          mask : pygac mask from reader
    Outputs:
      uncert : xarray dataset containing random and systematic uncertainty
               components
    """

    #
    # Get parameters for kernels/uncertainties
    #
    window, prt_bias, prt_sys, prt_threshold, ict_threshold, \
        space_threshold = get_uncert_parameter_thresholds()

    if ds["channels"].values.shape[1] == 409:
        gacdata = True
    else:
        gacdata = False

    avhrr_name = ds.attrs["spacecraft_name"]

    #
    # Get calibration coefficients and build channel specs
    #
    cal = Calibrator(ds.attrs["spacecraft_name"])
    specs = ir_channel_specs(cal, avhrr_name)
    twelve_micron = len(specs) == 3

    # Convenience aliases for the three possible channels (used throughout)
    spec_37, spec_11 = specs[0], specs[1]
    spec_12 = specs[2] if twelve_micron else None
    convT1 = spec_37.conv
    convT2 = spec_11.conv
    convT3 = spec_12.conv if twelve_micron else None
    NS_2, c0_2, c1_2, c2_2 = spec_11.space_radiance, *spec_11.nonlin_coeffs
    if twelve_micron:
        NS_3, c0_3, c1_3, c2_3 = spec_12.space_radiance, *spec_12.nonlin_coeffs

    #
    # Build per-channel data (noise + calibration variables) in one step
    #
    total_space = ds["full_space_counts"].values[:, :, :]
    total_ict = ds["full_ict_counts"].values[:, :, :]
    channels, bad_scan, Tict, ict1, ict2, ict3, ict4 = build_ir_channel_data(
        ds, specs, total_space, total_ict,
        window, prt_threshold, ict_threshold, space_threshold,
        gacdata, cal, mask,
    )

    # Convenience aliases kept for the rest of the function (pre-loop code,
    # per-scanline loop) — these will be removed as the loop is vectorised.
    ch_37, ch_11 = channels[0], channels[1]
    ch_12 = channels[2] if twelve_micron else None
    noise1, av_noise1, av_ict_noise1 = ch_37.noise, ch_37.av_noise, ch_37.av_ict_noise
    noise2, av_noise2, av_ict_noise2 = ch_11.noise, ch_11.av_noise, ch_11.av_ict_noise
    noise3, av_noise3, av_ict_noise3 = (
        (ch_12.noise, ch_12.av_noise, ch_12.av_ict_noise) if twelve_micron else (None, None, None)
    )
    CS_1, CICT_1, CE_1 = ch_37.cs, ch_37.cict, ch_37.ce
    CS_2, CICT_2, CE_2 = ch_11.cs, ch_11.cict, ch_11.ce
    if twelve_micron:
        CS_3, CICT_3, CE_3 = ch_12.cs, ch_12.cict, ch_12.ce

    solar_flag = np.zeros(CE_2.shape[0], dtype=np.uint8)
    if gacdata:
        #
        # See if solar contamination present
        # Only for GAC data
        # Possible at 2 locations
        #
        min_solar_1, max_solar_1, peak_solar_1, solar_solza_1,\
        min_solar_2, max_solar_2, peak_solar_2, solar_solza_2 = \
            find_solar(ds,mask,convT1)
        if min_solar_1 is None and max_solar_1 is None:
            min_solar_1 = -1
            max_solar_1 = -1
        if min_solar_2 is None and max_solar_2 is None:
            min_solar_2 = -1
            max_solar_2 = -1
    else:
        #
        # Find if stored solar contamination is present
        #
        min_solar_1, max_solar_1, min_solar_2, max_solar_2 = \
            get_solar_from_file(avhrr_name,ds)
    #
    # Set solar flag
    #
    if min_solar_1 >= 0 and max_solar_1 >= 0:
        solar_flag[min_solar_1:max_solar_1+1] = 1
    if min_solar_2 >= 0 and max_solar_2 >= 0:
        solar_flag[min_solar_2:max_solar_2+1] = 1

    #
    # Systematic components - uICT from gain variation in 3.7mu channel
    #
    gd = np.isfinite(ds["times"].values)
    time = ds["times"].values[gd][0]
    intimes = ds["times"].values[:]
    #
    # Only redo calculation if gac data
    #
    gain_37,gain_time = get_gainval(time,intimes,avhrr_name,ict1,ict2,
                                    ict3,ict4,CS_1,CICT_1,CE_1,0.,
                                    bad_scan,convT1,window,
                                    calculate=gacdata)
    if gain_37 is not None and gain_time is not None:
        uICT = get_uICT(gain_37,CS_1,CICT_1,Tict,0.,convT1,bad_scan,
                        solar_flag,window)
    else:
        #
        # Time of file out of gain time limits for HRPT data
        # Set uICT to NaNs
        #
        uICT = np.zeros(len(CS_1))
        uICT[:] = np.nan
    #
    # Vectorised computation — replaces the per-scanline for loop.
    # All operations are pure NumPy; (N,) arrays are reshaped to (N,1)
    # to broadcast against CE (N,P).
    #
    def col(x):
        """Reshape (N,) → (N,1) for broadcasting with CE (N,P)."""
        return np.asarray(x)[:, np.newaxis]

    ict_random1, ict_sys1, Lict_1 = get_ict_uncert(Tict, prt_bias, prt_sys, uICT, convT1)
    ict_random2, ict_sys2, Lict_2 = get_ict_uncert(Tict, prt_bias, prt_sys, uICT, convT2)
    if twelve_micron:
        ict_random3, ict_sys3, Lict_3 = get_ict_uncert(Tict, prt_bias, prt_sys, uICT, convT3)

    rad_37 = get_pixel(col(Lict_1), col(CS_1), CE_1, col(CICT_1), 0., 0., 0., 0.)
    rad_11 = get_pixel(col(Lict_2), col(CS_2), CE_2, col(CICT_2), NS_2, c0_2, c1_2, c2_2)
    if twelve_micron:
        rad_12 = get_pixel(col(Lict_3), col(CS_3), CE_3, col(CICT_3), NS_3, c0_3, c1_3, c2_3)

    rad_noise_37 = _radiance_random_uncert(
        spec_37, noise1, av_noise1, av_ict_noise1,
        col(ict_random1), col(Lict_1), col(CS_1), CE_1, col(CICT_1))
    rad_noise_11 = _radiance_random_uncert(
        spec_11, noise2, av_noise2, av_ict_noise2,
        col(ict_random2), col(Lict_2), col(CS_2), CE_2, col(CICT_2))
    if twelve_micron:
        rad_noise_12 = _radiance_random_uncert(
            spec_12, noise3, av_noise3, av_ict_noise3,
            col(ict_random3), col(Lict_3), col(CS_3), CE_3, col(CICT_3))

    _, bt_rand_37 = convT1.rad_to_t_uncert(rad_37, rad_noise_37)
    _, bt_rand_11 = convT2.rad_to_t_uncert(rad_11, rad_noise_11)
    if twelve_micron:
        _, bt_rand_12 = convT3.rad_to_t_uncert(rad_12, rad_noise_12)
    else:
        bt_rand_12 = np.full(CE_2.shape, np.nan, dtype=CE_2.dtype)

    # Systematic: NaN propagates naturally where ict_sys / uICT is NaN.
    rad_sys_37 = _radiance_sys_uncert(spec_37, col(ict_sys1), col(Tict), col(CS_1), CE_1, col(CICT_1))
    rad_sys_11 = _radiance_sys_uncert(spec_11, col(ict_sys2), col(Tict), col(CS_2), CE_2, col(CICT_2))
    if twelve_micron:
        rad_sys_12 = _radiance_sys_uncert(spec_12, col(ict_sys3), col(Tict), col(CS_3), CE_3, col(CICT_3))

    _, bt_sys_37 = convT1.rad_to_t_uncert(rad_37, rad_sys_37)
    _, bt_sys_11 = convT2.rad_to_t_uncert(rad_11, rad_sys_11)
    if twelve_micron:
        _, bt_sys_12 = convT3.rad_to_t_uncert(rad_12, rad_sys_12)
    else:
        bt_sys_12 = np.full(CE_2.shape, np.nan, dtype=CE_2.dtype)

    # Add measurement equation uncertainty (0.5K@300K for 3.7µm; flat 0.5K for 11/12µm).
    delta_rad_37 = (convT1.t_to_rad(300. + 0.5 / np.sqrt(3.))
                    - convT1.t_to_rad(300. - 0.5 / np.sqrt(3.))) / 2.
    _, new_37_uncert = convT1.rad_to_t_uncert(rad_37, delta_rad_37)
    tot_sys_37 = np.sqrt(bt_sys_37 ** 2 + new_37_uncert ** 2)
    tot_sys_11 = np.sqrt(bt_sys_11 ** 2 + 0.5 ** 2 / 3.)
    if twelve_micron:
        tot_sys_12 = np.sqrt(bt_sys_12 ** 2 + 0.5 ** 2 / 3.)
    else:
        tot_sys_12 = np.full(CE_2.shape, np.nan, dtype=CE_2.dtype)

    uratio_37 = bt_sys_37 / tot_sys_37   # NaN/NaN = NaN where no valid sys ✓
    uratio_11 = bt_sys_11 / tot_sys_11
    if twelve_micron:
        uratio_12 = bt_sys_12 / tot_sys_12
    else:
        uratio_12 = np.full(CE_2.shape, np.nan, dtype=CE_2.dtype)

    bt_sys_37 = tot_sys_37
    bt_sys_11 = tot_sys_11
    bt_sys_12 = tot_sys_12

    # Apply bad-scan mask last (equivalent to the `continue` branch in the old loop).
    bad = bad_scan == 1
    for arr in (bt_rand_37, bt_rand_11, bt_rand_12,
                bt_sys_37, bt_sys_11, bt_sys_12,
                uratio_37, uratio_11, uratio_12):
        arr[bad] = np.nan

    #
    # Output uncertainties
    #
    random = np.zeros((bt_rand_11.shape[0],bt_rand_11.shape[1],3))
    systematic = np.zeros((bt_rand_11.shape[0],bt_rand_11.shape[1],3))
    uratio = np.zeros((bt_rand_11.shape[0],bt_rand_11.shape[1],3),
                      dtype=np.uint8)
    uflags = np.zeros((bt_rand_11.shape[0]),dtype=np.uint8)

    random[:,:,0] = bt_rand_37
    random[:,:,1] = bt_rand_11
    if twelve_micron:
        random[:,:,2] = bt_rand_12
    else:
        random[:,:,2] = np.nan
    systematic[:,:,0] = bt_sys_37
    systematic[:,:,1] = bt_sys_11
    if twelve_micron:
        systematic[:,:,2] = bt_sys_12
    else:
        systematic[:,:,2] = np.nan
    #
    # Ratio for channel-to-channel covariance as ubyte
    #
    uratio[:,:,:] = 0.
    gd = np.isfinite(uratio_37)&(uratio_37 < 0.)
    uratio_37[gd] = 0.
    gd = np.isfinite(uratio_37)&(uratio_37 > 1.)
    uratio_37[gd] = 1.
    gd = np.isfinite(uratio_11)&(uratio_11 < 0.)
    uratio_11[gd] = 0.
    gd = np.isfinite(uratio_11)&(uratio_11 > 1.)
    uratio_11[gd] = 1.
    gd = np.isfinite(uratio_37)
    uratio[gd,0] = (uratio_37[gd]*255).astype(dtype=np.uint8)
    gd = np.isfinite(uratio_11)
    uratio[gd,1] = (uratio_11[gd]*255).astype(dtype=np.uint8)
    if twelve_micron:
        gd = np.isfinite(uratio_12)&(uratio_12 < 0.)
        uratio_12[gd] = 0.
        gd = np.isfinite(uratio_12)&(uratio_12 > 1.)
        uratio_12[gd] = 1.
        gd = np.isfinite(uratio_12)
        uratio[gd,2] = (uratio_12[gd]*255).astype(dtype=np.uint8)
    else:
        uratio[:,:,2] = 0

    #
    # Flags
    #
    gd = (bad_scan == 1)
    uflags[gd] = 1
    gd = (solar_flag == 1)
    uflags[gd] = (uflags[gd]|2)
    if np.sum(np.isfinite(systematic[:,:,1])) == 0:
        uflags = (uflags|4)

    time = (ds["times"].values - np.datetime64("1970-01-01 00:00:00"))/\
           np.timedelta64(1,"s")
    time_da = xr.DataArray(time,dims=["times"],attrs={"long_name":"scanline time",
                                                     "units":"seconds since 1970-01-01"})
    across_da = xr.DataArray(np.arange(random.shape[1]),dims=["across_track"])
    ir_channels_da = xr.DataArray(np.array([3,4,5]),dims=["ir_channels"])
    random_da = xr.DataArray(random,
                             dims=["times","across_track","ir_channels"],
                             attrs={"long_name":"Random uncertainties","units":"K"})
    sys_da = xr.DataArray(systematic,
                          dims=["times","across_track","ir_channels"],
                          attrs={"long_name":"Systematic uncertainties","units":"K"})

    uratio_da = xr.DataArray(uratio,
                             dims=["times","across_track","ir_channels"],
                             attrs={"long_name":"Channel-to-channel covariance  ratio",
                                    "_FillValue":0})

    uflags_da = xr.DataArray(uflags,
                             dims=["times"],
                             attrs={"long_name":"Uncertainty flags",
                                    "flag_masks": "1b, 2b, 4b",
                                    "flag_meanings": ("bad_space_view "
                                                      "solar_contamination_of_gain "
                                                      "no_IR_systematic_uncertainty ")})

    uncertainties = xr.Dataset(dict(times=time_da,across_track=across_da,
                                    ir_channels=ir_channels_da,
                                    random=random_da,systematic=sys_da,
                                    chan_covar_ratio=uratio_da,
                                    uncert_flags=uflags_da))

    return uncertainties
