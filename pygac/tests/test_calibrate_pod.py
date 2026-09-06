#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014, 2015 Martin Raspaud

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>

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

"""Test function for the POD calibration.
"""


import unittest

import numpy as np
import pytest

from pygac.calibration.noaa import Calibrator, calibrate_solar, calibrate_thermal


class TestGenericCalibration(unittest.TestCase):

    def test_calibration_vis(self):

        counts = np.array([[0, 0, 0, 0, 0,
                            512, 512, 512, 512, 512,
                            1023, 1023, 1023, 1023, 1023],
                           [41, 41, 41, 41, 41,
                            150, 150, 150, 150, 150,
                            700, 700, 700, 700, 700]])
        year = 1997
        jday = 196
        spacecraft_id = "noaa14"
        cal = Calibrator(spacecraft_id)
        corr = 1

        channel = 0

        ref1 = calibrate_solar(counts[:, channel::5], channel, year, jday, cal, corr)

        channel = 1

        ref2 = calibrate_solar(counts[:, channel::5], channel, year, jday, cal, corr)

        channel = 2

        data = np.ma.array(counts[:, channel::5], mask=True)

        ref3 = calibrate_solar(data, channel, year, jday, cal, corr)

        expected = (np.array([[np.nan, 60.91525491, 127.00377987],
                              [0., 14.0971609, 85.22962417]]),
                    np.array([[np.nan, 72.51635437, 151.19121018],
                              [0., 16.7819164, 101.46131111]]),
                    np.array([[-32001., -32001., -32001.],
                              [-32001., -32001., -32001.]]))

        np.testing.assert_allclose(ref1, expected[0])
        np.testing.assert_allclose(ref2, expected[1])
        np.testing.assert_allclose(ref3.filled(-32001), expected[2])

    def test_calibration_ir(self):
        counts = np.array([[0, 0, 612, 0, 0,
                            512, 512, 487, 512, 512,
                            923, 923, 687, 923, 923],
                           [41, 41, 634, 41, 41,
                            150, 150, 461, 150, 150,
                            700, 700, 670, 700, 700],
                           [241, 241, 656, 241, 241,
                            350, 350, 490, 350, 350,
                            600, 600, 475, 600, 600]])
        prt_counts = np.array([0, 230, 230])
        ict_counts = np.array([[745.3, 397.9, 377.8],
                               [744.8, 398.1, 378.4],
                               [745.7, 398., 378.3]])
        space_counts = np.array([[987.3, 992.5, 989.4],
                                 [986.9, 992.8, 989.6],
                                 [986.3, 992.3, 988.9]])

        spacecraft_id = "noaa14"
        cal = Calibrator(spacecraft_id)
        ch3 = calibrate_thermal(counts[:, 2::5],
                                prt_counts,
                                ict_counts[:, 0],
                                space_counts[:, 0],
                                line_numbers=np.array([1, 2, 3]),
                                channel=3,
                                cal=cal)

        expected_ch3 = np.array([[298.28466, 305.167571, 293.16182],
                                 [296.878502, 306.414234, 294.410224],
                                 [295.396779, 305.020259, 305.749526]])

        np.testing.assert_allclose(expected_ch3, ch3)

        ch4 = calibrate_thermal(counts[:, 3::5],
                                prt_counts,
                                ict_counts[:, 1],
                                space_counts[:, 1],
                                line_numbers=np.array([1, 2, 3]),
                                channel=4,
                                cal=cal)

        expected_ch4 = np.array([[325.828062, 275.414804, 196.214709],
                                 [322.359517, 312.785057, 249.380649],
                                 [304.326806, 293.490822, 264.148021]])

        np.testing.assert_allclose(expected_ch4, ch4)

        ch5 = calibrate_thermal(counts[:, 4::5],
                                prt_counts,
                                ict_counts[:, 2],
                                space_counts[:, 2],
                                line_numbers=np.array([1, 2, 3]),
                                channel=5,
                                cal=cal)

        expected_ch5 = np.array([[326.460316, 272.146547, 187.434456],
                                 [322.717606, 312.388155, 244.241633],
                                 [303.267012, 291.590832, 260.05426]])

        np.testing.assert_allclose(expected_ch5, ch5)


class TestThermalCalibrationOnRealData:
    """Pin the thermal calibration output on a real POD scene.

    ``test_calibration_ir`` above exercises only three scan lines, which is
    below the 51-line smoothing window and so never reaches the convolution
    and broadcasting path used for real passes. This characterisation test
    covers it, so that refactoring ``calibrate_thermal`` cannot silently move
    brightness temperatures.
    """

    # channels 3, 4 and 5 at [::97, ::41], from the committed fixture.
    # None marks pixels the fixture masks out.
    EXPECTED_IR = [
        [[313.7021087861816, 301.3761485641495, 300.5297275608002],
         [314.40218317509635, 304.0959679404656, 303.4547465036669],
         [314.9170859410342, 304.8036868338647, 304.1817618772081],
         [314.9798074251143, 305.51290615962046, 304.98203016360105]],
        [[313.9380064172301, 302.8260896860293, 302.1393818282739],
         [314.7131517129142, 304.78336011232403, 304.2504874227244],
         [None, 305.92848505218126, 305.4437994203835],
         [314.2157300965453, 304.122312584772, 303.4768212783679]],
        [[None, None, None], [None, None, None],
         [None, None, None], [None, None, None]],
        [[314.819002923425, 305.09913870821595, 304.52452744874347],
         [314.8811341295552, 305.93046036091704, 305.44524425617413],
         [315.12874234743674, 306.32163098655764, 305.824829264994],
         [314.9513754511419, 305.9540839782427, 305.4946167227422]],
    ]

    def test_thermal_calibration_is_unchanged(self):
        from pathlib import Path

        import xarray as xr

        from pygac.calibration.noaa import calibrate

        path = (Path(__file__).parent / "data" / "uncertainty_regression"
                / "noaa14_pod_d00322.input.nc")
        if not path.exists():
            pytest.skip("quick-tier fixture not found")
        with xr.open_dataset(path) as ds:
            channels = calibrate(ds.load())["channels"].values

        got = channels[::97, ::41, 2:]
        expected = np.array(
            [[[np.nan if v is None else v for v in col] for col in row]
             for row in self.EXPECTED_IR]
        )
        np.testing.assert_allclose(got, expected, rtol=1e-12, equal_nan=True)
