#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014-2019 Pytroll Developers

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
try:
    import mock
except ImportError:
    from unittest import mock
import numpy as np

from pygac.calibration import Calibrator, calibrate_solar, calibrate_thermal


class TestGenericCalibration(unittest.TestCase):

    def test_calibration_vis(self):

        counts = np.array([[0, 0, 0, 0, 0,
                            512, 512, 512, 512, 512,
                            1023, 1023, 1023, 1023, 1023],
                           [41, 41, 41, 41, 41,
                            150, 150, 150, 150, 150,
                            700, 700, 700, 700, 700]])
        year = 2010
        jday = 1
        spacecraft_id = "noaa19"
        cal = Calibrator(spacecraft_id)
        corr = 1
        channel = 0

        ref1 = calibrate_solar(counts[:, channel::5], channel, year, jday, cal, corr)

        channel = 1

        ref2 = calibrate_solar(counts[:, channel::5], channel, year, jday, cal, corr)

        channel = 2

        data = np.ma.array(counts[:, channel::5], mask=True)

        ref3 = calibrate_solar(data, channel, year, jday, cal, corr)

        expected = (np.array([[np.nan, 27.37909518, 110.60103456],
                              [0.11943135, 6.03671211, 57.99695154]]),
                    np.array([[np.nan, 3.05229160e+01, 1.24811455e+02],
                              [1.23011792e-01, 6.82715447e+00, 6.52122414e+01]]),
                    np.array([[0., 523.41775, 1034.41775],
                              [41., 150., 711.41775]]))
        np.testing.assert_allclose(ref1, expected[0])
        np.testing.assert_allclose(ref2, expected[1])
        np.testing.assert_allclose(ref3, expected[2])

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
        space_counts = np.array([[987.3,  992.5,  989.4],
                                 [986.9,  992.8,  989.6],
                                 [986.3,  992.3,  988.9]])

        spacecraft_id = "noaa19"
        cal = Calibrator(spacecraft_id)
        ch3 = calibrate_thermal(counts[:, 2::5],
                                prt_counts,
                                ict_counts[:, 0],
                                space_counts[:, 0],
                                line_numbers=np.array([1, 2, 3]),
                                channel=3,
                                cal=cal)

        expected_ch3 = np.array([[298.36742, 305.248478, 293.238328],
                                 [296.960275, 306.493766, 294.488956],
                                 [295.476935, 305.101309, 305.829827]])

        np.testing.assert_allclose(expected_ch3, ch3)

        ch4 = calibrate_thermal(counts[:, 3::5],
                                prt_counts,
                                ict_counts[:, 1],
                                space_counts[:, 1],
                                line_numbers=np.array([1, 2, 3]),
                                channel=4,
                                cal=cal)

        expected_ch4 = np.array([[326.576534, 275.348988, 197.688755],
                                 [323.013104, 313.207077, 249.36352],
                                 [304.58091, 293.579308, 264.0631]])

        np.testing.assert_allclose(expected_ch4, ch4)

        ch5 = calibrate_thermal(counts[:, 4::5],
                                prt_counts,
                                ict_counts[:, 2],
                                space_counts[:, 2],
                                line_numbers=np.array([1, 2, 3]),
                                channel=5,
                                cal=cal)

        expected_ch5 = np.array([[326.96161, 272.090164, 188.267991],
                                 [323.156317, 312.673269, 244.184452],
                                 [303.439383, 291.649444, 259.973091]])

        np.testing.assert_allclose(expected_ch5, ch5)
