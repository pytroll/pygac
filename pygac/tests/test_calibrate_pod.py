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

from pygac.calibration import calibrate_solar, calibrate_thermal


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
        channel3_switch = 0
        corr = 1
        number_of_data_records = 2

        channel = 0

        ref1 = calibrate_solar(counts[:, channel::5], channel, year, jday,
                               spacecraft_id, corr)

        channel = 1

        ref2 = calibrate_solar(counts[:, channel::5], channel, year, jday,
                               spacecraft_id, corr)

        channel = 2

        data = np.ma.array(counts[:, channel::5], mask=True)

        ref3 = calibrate_solar(data, channel, year, jday,
                               spacecraft_id, corr)

        expected = (np.array([[-5.34517394,   61.40431529,  128.02343443],
                              [0.,   14.21034048,   85.91389337]]),
                    np.array([[-6.50301792,   74.70540093,  155.75520959],
                              [0.,   17.28851104,  104.52411723]]),
                    np.array([[-32001., -32001., -32001.],
                              [-32001., -32001., -32001.]]))

        self.assertTrue(np.allclose(ref1, expected[0]))
        self.assertTrue(np.allclose(ref2, expected[1]))
        self.assertTrue(np.allclose(ref3.filled(-32001), expected[2]))

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

        spacecraft_id = "noaa14"
        number_of_data_records = 3
        ch3 = calibrate_thermal(counts[:, 2::5],
                                prt_counts,
                                ict_counts[:, 0],
                                space_counts[:, 0],
                                line_numbers=np.array([1, 2, 3]),
                                channel=3,
                                spacecraft=spacecraft_id)

        expected_ch3 = np.array([[298.28524223, 305.16852862, 293.16212655],
                                 [296.87900835, 306.41526012, 294.41059746],
                                 [295.39720547, 305.02120845, 305.75051609]])

        self.assertTrue(np.allclose(expected_ch3, ch3))

        ch4 = calibrate_thermal(counts[:, 3::5],
                                prt_counts,
                                ict_counts[:, 1],
                                space_counts[:, 1],
                                line_numbers=np.array([1, 2, 3]),
                                channel=4,
                                spacecraft=spacecraft_id)

        expected_ch4 = np.array([[325.82572316, 275.41391816, 196.21443457],
                                 [322.35731456, 312.78320066, 249.38013251],
                                 [304.32522137, 293.48953883, 264.14732182]])

        self.assertTrue(np.allclose(expected_ch4, ch4))

        ch5 = calibrate_thermal(counts[:, 4::5],
                                prt_counts,
                                ict_counts[:, 2],
                                space_counts[:, 2],
                                line_numbers=np.array([1, 2, 3]),
                                channel=5,
                                spacecraft=spacecraft_id)

        expected_ch5 = np.array([[326.47287181,  272.14169523,  187.40907142],
                                 [322.72885806,  312.39588991,  244.22910864],
                                 [303.27173737,  291.59183911,  260.0459766]])

        self.assertTrue(np.allclose(expected_ch5, ch5))


def suite():
    """The suite for test_slerp
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestGenericCalibration))

    return mysuite


if __name__ == '__main__':
    unittest.main()
