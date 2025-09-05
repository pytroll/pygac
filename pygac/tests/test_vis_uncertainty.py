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
import xarray as xr
import pandas as pd
from pygac.calibration.ir_uncertainty import allan_deviation, get_bad_space_counts, get_uncert_parameter_thresholds
from pygac.calibration.noaa import Calibrator
from pygac.calibration.vis_uncertainty import get_vars, \
                                                get_gain, get_random, get_sys, get_noise, vis_uncertainty

n_scan_lines = 100
n_columns = 1
n_channels = 5
n_ir_channels = 3
n_vis_channels = 3
n_pixels = 10

# Coordinate values
scan_line_index = np.arange(1, n_scan_lines + 1, dtype=np.int16)
columns = np.arange(n_columns, dtype=np.int32)
channel_name = np.array(['1', '2', '3', '4', '5'], dtype='<U1')
ir_channel_name = np.array(['3', '4', '5'], dtype='<U1')
vis_channel_name = np.array(['1', '2', '3'], dtype='<U1')
pixel_index = np.arange(n_pixels, dtype=np.int8)
times = pd.date_range("1987-02-02", periods=n_scan_lines, freq="S")

# Sample data
longitude = np.random.uniform(-180, 180, (n_scan_lines, n_columns)).astype(np.float32)
latitude = np.random.uniform(-90, 90, (n_scan_lines, n_columns)).astype(np.float32)
channels = np.random.rand(n_scan_lines, n_columns, n_channels)
counts = np.random.rand(n_scan_lines, n_pixels, n_vis_channels)
vis_space_counts = np.random.rand(n_scan_lines, n_vis_channels)
total_vis_space_counts = np.random.rand(n_scan_lines, n_pixels, n_vis_channels)
sun_zen = np.random.uniform(0, 180, (n_scan_lines, n_columns)).astype(np.float32)


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
        "counts": (["scan_line_index", "pixel_index", "vis_channel_name"], counts),
        "vis_space_counts": (["scan_line_index", "vis_channel_name"], vis_space_counts),
        "total_vis_space_counts": (["scan_line_index", "pixel_index", "vis_channel_name"], total_vis_space_counts),
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

    # def test_FOV_solar_contam(self):
    #     refl = np.array([0.025, 0.055, 0.07, 0.03, 0.067])
    #     sza = np.array([101, 106, 99, 104, 105])
    #
    #     window, solar_contam_threshold, sza_threshold = \
    #         get_uncert_parameter_thresholds(vischans=True)
    #
    #     fov_meas = get_FOV_solar_contam(refl, sza, 0.05, 102.)
    #     self.assertEqual(np.count_nonzero(fov_meas), 2)

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
