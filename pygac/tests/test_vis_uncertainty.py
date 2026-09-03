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

"""Unit tests for uncertainty handling in the VIS channels.
"""

import unittest

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from pygac.calibration.noaa import Calibrator
from pygac.uncertainty.ir import allan_deviation, get_bad_space_counts
from pygac.uncertainty.vis import get_gain, get_sys, get_vars

n_scan_lines = 100
n_columns = 1
n_channels = 5
n_ir_channels = 3
n_vis_channels = 3
n_pixels = 10

# Coordinate values
scan_line_index = np.arange(1, n_scan_lines + 1, dtype=np.int16)
columns = np.arange(n_columns, dtype=np.int32)
channel_name = np.array(["1", "2", "3", "4", "5"], dtype="<U1")
ir_channel_name = np.array(["3", "4", "5"], dtype="<U1")
vis_channel_name = np.array(["1", "2", "3"], dtype="<U1")
pixel_index = np.arange(n_pixels, dtype=np.int8)
times = pd.date_range("1987-02-02", periods=n_scan_lines, freq="S")

# Sample data
rng = np.random.default_rng()

longitude = rng.uniform(-180, 180, (n_scan_lines, n_columns)).astype(np.float32)
latitude = rng.uniform(-90, 90, (n_scan_lines, n_columns)).astype(np.float32)
channels = rng.random((n_scan_lines, n_columns, n_channels))
counts = rng.random((n_scan_lines, n_columns, n_channels))
full_space_counts = rng.random((n_scan_lines, n_pixels, n_channels))
sun_zen = rng.uniform(0, 180, (n_scan_lines, n_columns)).astype(np.float32)


# Construct the Dataset
ds = xr.Dataset(
    coords={
        "scan_line_index": ("scan_line_index", scan_line_index),
        "columns": ("columns", columns),
        "channel_name": ("channel_name", channel_name),
        "vis_channel_name": ("vis_channel_name", vis_channel_name),
        "pixel_index": ("pixel_index", pixel_index),
        "times": ("scan_line_index", times),
        "longitude": (["scan_line_index", "columns"], longitude),
        "latitude": (["scan_line_index", "columns"], latitude),
    },
    data_vars={
        "channels": (["scan_line_index", "columns", "channel_name"], channels),
        "counts": (["scan_line_index", "columns", "channel_name"], counts),
        "full_space_counts": (["scan_line_index", "pixel_index", "channel_name"], full_space_counts),
        "sun_zen": (["scan_line_index", "columns"], sun_zen),
        },
    attrs={
        "spacecraft_name": "noaa10",
        "sun_earth_distance_correction_factor": 0.9712376984147952,
    }
)


class TestVisibleUncertainty(unittest.TestCase):
    def test_Allan_deviation(self):
        self.space = np.array([[37., 37., 36., 37., 38., 37., 37., 36., 38., 37.],
                          [36., 37., 35., 37., 38., 38., 37., 36., 36., 37.]])
        measurement = allan_deviation(self.space, bad_scan=None)

        exp_measurement = 0.79930525

        self.assertAlmostEqual(measurement, exp_measurement)

    def test_bad_space_counts(self):
        sp_data = np.array([[37., 37., 36., 37., 42., 37., 37., 36., 28., 37.],
                          [36., 37., 25., 37., 38., 45., 37., 36., 36., 37.]])

        bad_counts = get_bad_space_counts(sp_data)

        self.assertEqual(np.count_nonzero(bad_counts), 4)

    def test_vis_uncert(self):
        cal = Calibrator(ds.attrs["spacecraft_name"])
        s0_1 = cal.s0[0]
        s1_1 = cal.s1[0]
        s2_1 = cal.s2[0]
        times = ds.coords["times"]
        start_time = times[0].dt
        year = start_time.year.item()
        jday = start_time.dayofyear.item()
        l_date = Calibrator.date2float(cal.date_of_launch)
        time = (year + jday / 365.0) - l_date
        gain_1 = get_gain(s0_1, s1_1, s2_1, time, cal, 0)
        #
        # Check gain
        #
        check_gain = 0.11336
        np.testing.assert_allclose(gain_1, check_gain, atol=0.00001)

        D_1, C_1 = get_vars(ds, 0)
        u_sys_exp = np.zeros(C_1.shape, dtype=C_1.dtype)
        for i in range(len(D_1)):
            u_sys_exp[i,:] = ((C_1[i,:] - D_1[i])**2)*((0.050990195135927854*gain_1)**2)
            u_sys_exp[i,:] = np.sqrt(u_sys_exp[i,:])

        u_sys = np.zeros(C_1.shape, dtype=C_1.dtype)
        for i in range(len(D_1)):
            u_sys[i,:] = get_sys(1, C_1[i,:], D_1[i], gain_1)

        np.testing.assert_allclose(u_sys, u_sys_exp, atol=0.0001)


class TestComputeVisBadScans:
    """_compute_vis_bad_scans vectorises the bad-scan detection loop."""

    def _space_all_good(self, n=5, p=10, nc=2):
        """All values near mean — no bad pixels."""
        rng = np.random.default_rng(0)
        return (37.0 + rng.uniform(-0.5, 0.5, (n, p, nc))).astype(float)

    def test_no_bad_scans_when_all_good(self):
        from pygac.uncertainty.vis import _compute_vis_bad_scans
        space = self._space_all_good()
        bad = _compute_vis_bad_scans(space, chan_3a=False)
        assert np.all(bad == 0)

    def test_returns_int8(self):
        from pygac.uncertainty.vis import _compute_vis_bad_scans
        space = self._space_all_good()
        bad = _compute_vis_bad_scans(space, chan_3a=False)
        assert bad.dtype == np.int8

    def test_scanline_flagged_when_channel1_bad(self):
        from pygac.uncertainty.vis import _compute_vis_bad_scans
        space = self._space_all_good(n=4, p=10, nc=2)
        space[2, 3, 0] = 999.0   # outlier in ch0 scanline 2
        bad = _compute_vis_bad_scans(space, chan_3a=False)
        assert bad[2] == 1
        assert bad[0] == 0 and bad[1] == 0 and bad[3] == 0

    def test_chan3a_flagged_when_third_channel_bad(self):
        from pygac.uncertainty.vis import _compute_vis_bad_scans
        space = self._space_all_good(n=4, p=10, nc=3)
        space[1, 5, 2] = 999.0
        bad = _compute_vis_bad_scans(space, chan_3a=True)
        assert bad[1] == 1
        assert bad[0] == 0


class TestVisChannelNoise:
    """_vis_channel_noise returns (noise, av_noise) for a single channel."""

    def _good_space(self, n=10, p=10):
        rng = np.random.default_rng(1)
        return (37.0 + rng.uniform(-0.5, 0.5, (n, p))).astype(float)

    def test_returns_two_scalars(self):
        from pygac.uncertainty.vis import _vis_channel_noise
        space = self._good_space()
        bad = np.zeros(10, dtype=np.int8)
        noise, av_noise = _vis_channel_noise(space, bad, window=3)
        assert np.ndim(noise) == 0
        assert np.ndim(av_noise) == 0

    def test_av_noise_smaller_than_noise(self):
        from pygac.uncertainty.vis import _vis_channel_noise
        space = self._good_space()
        bad = np.zeros(10, dtype=np.int8)
        noise, av_noise = _vis_channel_noise(space, bad, window=3)
        assert av_noise < noise

    def test_digitisation_contribution(self):
        """noise must be >= sqrt(1/3) even for zero Allan deviation."""
        from pygac.uncertainty.vis import _vis_channel_noise
        space = np.full((10, 10), 37.0)   # flat → Allan dev = 0
        bad = np.zeros(10, dtype=np.int8)
        noise, _ = _vis_channel_noise(space, bad, window=3)
        assert noise >= np.sqrt(1.0 / 3)


class TestVisSysUncert:
    """_vis_sys_uncert replaces get_sys(channel, ...) with explicit bool flag."""

    def _args(self):
        rng = np.random.default_rng(2)
        counts = rng.uniform(200, 300, 128)
        mean_space = 37.0
        gain = 0.113
        return counts, mean_space, gain

    def test_without_wv_matches_get_sys_ch1(self):
        from pygac.uncertainty.vis import _vis_sys_uncert, get_sys
        counts, mean_space, gain = self._args()
        expected = get_sys(1, counts, mean_space, gain)
        result = _vis_sys_uncert(counts, mean_space, gain, include_water_vapour=False)
        np.testing.assert_allclose(result, expected)

    def test_with_wv_matches_get_sys_ch2(self):
        from pygac.uncertainty.vis import _vis_sys_uncert, get_sys
        counts, mean_space, gain = self._args()
        expected = get_sys(2, counts, mean_space, gain)
        result = _vis_sys_uncert(counts, mean_space, gain, include_water_vapour=True)
        np.testing.assert_allclose(result, expected)

    def test_with_wv_larger_than_without(self):
        from pygac.uncertainty.vis import _vis_sys_uncert
        counts, mean_space, gain = self._args()
        no_wv = _vis_sys_uncert(counts, mean_space, gain, include_water_vapour=False)
        wv = _vis_sys_uncert(counts, mean_space, gain, include_water_vapour=True)
        assert np.all(wv >= no_wv)


class TestVisRandomUncert:
    """_vis_random_uncert returns (uncert, Rcal) without unused cal/year/jday args."""

    def _args(self):
        rng = np.random.default_rng(3)
        noise, av_noise = 1.2, 0.38
        gain = 0.113
        counts = rng.uniform(200, 300, 128)
        mean_space = 37.0
        return noise, av_noise, gain, counts, mean_space

    def test_matches_get_random(self):
        from pygac.calibration.noaa import Calibrator
        from pygac.uncertainty.vis import _vis_random_uncert, get_random
        noise, av_noise, gain, counts, mean_space = self._args()
        cal = Calibrator("noaa10")
        exp_uncert, exp_rcal = get_random(noise, av_noise, gain, cal, 1987, 33, counts, mean_space)
        uncert, rcal = _vis_random_uncert(noise, av_noise, gain, counts, mean_space)
        np.testing.assert_allclose(uncert, exp_uncert)
        np.testing.assert_allclose(rcal, exp_rcal)

    def test_rcal_shape(self):
        from pygac.uncertainty.vis import _vis_random_uncert
        noise, av_noise, gain, counts, mean_space = self._args()
        _, rcal = _vis_random_uncert(noise, av_noise, gain, counts, mean_space)
        assert rcal.shape == counts.shape


class TestVisUncertaintyIntegration:
    """Behavioral tests for vis_uncertainty output (shape, masking, dtypes)."""

    @pytest.fixture
    def synthetic_ds(self):
        import pandas as pd
        import xarray as xr
        rng = np.random.default_rng(42)
        n, p, nc = 20, 128, 5
        sp = np.zeros((n, 10, nc))
        sp[:, :, :] = 37.0 + rng.uniform(-0.3, 0.3, (n, 10, nc))
        # Inject a bad space pixel in scanline 5 channel 0
        sp[5, 3, 2] = 999.0
        times = pd.date_range("2000-01-01", periods=n, freq="s")
        return xr.Dataset(
            data_vars={
                "channels": (["times", "columns", "ch"], rng.random((n, p, nc))),
                "counts": (["times", "columns", "vis"], rng.random((n, p, 2))),
                "full_space_counts": (["times", "pixels", "ch"], sp),
                "sun_zen": (["times", "columns"], rng.uniform(0, 80, (n, p)).astype(np.float32)),
            },
            coords={"times": times},
            attrs={
                "spacecraft_name": "noaa14",
                "sun_earth_distance_correction_factor": 0.97,
            },
        )

    def test_output_is_dataset(self, synthetic_ds):
        import xarray as xr

        from pygac.uncertainty.vis import vis_uncertainty
        result = vis_uncertainty(synthetic_ds, mask=None)
        assert isinstance(result, xr.Dataset)

    def test_random_shape(self, synthetic_ds):
        from pygac.uncertainty.vis import vis_uncertainty
        result = vis_uncertainty(synthetic_ds, mask=None)
        assert result["random"].shape == (20, 128, 3)

    def test_bad_scan_produces_nan_random(self, synthetic_ds):
        from pygac.uncertainty.vis import vis_uncertainty
        result = vis_uncertainty(synthetic_ds, mask=None)
        # scanline 5 should be NaN (bad space view)
        assert np.all(np.isnan(result["random"].values[5, :, 0]))
        # scanline 0 should be finite
        assert np.all(np.isfinite(result["random"].values[0, :, 0]))
