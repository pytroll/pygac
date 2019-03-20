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
        year = 2010
        jday = 1
        spacecraft_id = "noaa19"
        corr = 1
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

        expected = (np.array([[-2.11703668, 27.53393535, 111.69503509],
                              [0.12003816, 6.06738347, 58.49731451]]),
                    np.array([[-2.40976061e+00, 3.06514966e+01, 1.24856018e+02],
                              [1.23577467e-01, 6.85854943e+00, 6.53099114e+01]]),
                    np.array([[0.,   512.,  1023.],
                              [41.,   150.,   700.]]))
        self.assertTrue(np.allclose(ref1, expected[0]))
        self.assertTrue(np.allclose(ref2, expected[1]))
        self.assertTrue(np.allclose(ref3, expected[2]))

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
        ch3 = calibrate_thermal(counts[:, 2::5],
                                prt_counts,
                                ict_counts[:, 0],
                                space_counts[:, 0],
                                line_numbers=np.array([1, 2, 3]),
                                channel=3,
                                spacecraft=spacecraft_id)

        expected_ch3 = np.array([[298.36772477, 305.24899954, 293.23847375],
                                 [296.96053595, 306.49432811, 294.48914038],
                                 [295.47715016, 305.10182601, 305.83036782]])

        self.assertTrue(np.allclose(expected_ch3, ch3))

        ch4 = calibrate_thermal(counts[:, 3::5],
                                prt_counts,
                                ict_counts[:, 1],
                                space_counts[:, 1],
                                line_numbers=np.array([1, 2, 3]),
                                channel=4,
                                spacecraft=spacecraft_id)

        expected_ch4 = np.array([[326.57669548, 275.34893211, 197.68844955],
                                 [323.01324859, 313.20717645, 249.3633716],
                                 [304.58097221, 293.57932356, 264.0630027]])

        self.assertTrue(np.allclose(expected_ch4, ch4))

        ch5 = calibrate_thermal(counts[:, 4::5],
                                prt_counts,
                                ict_counts[:, 2],
                                space_counts[:, 2],
                                line_numbers=np.array([1, 2, 3]),
                                channel=5,
                                spacecraft=spacecraft_id)

        expected_ch5 = np.array([[326.96168274, 272.09013413, 188.26784127],
                                 [323.15638147, 312.67331324, 244.18437795],
                                 [303.43940924, 291.64944851, 259.97304154]])

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
