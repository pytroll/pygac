#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014-2019 Pytroll Developers

# Author(s):

#   Nina Hakansson <nina.hakansson@smhi.se>
#   Adam Dybbroe <adam.dypbroe@smhi.se>

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

"""Test function for the angle calculation.
"""

import unittest

import numpy as np

from pygac import get_absolute_azimuth_angle_diff


class TestAngles(unittest.TestCase):

    def test_azidiff_angle(self):


        SAT_AZ = np.ma.array([[48.0,  56.0, 64.0,  72.0],
                              [80.0,  88.0, 96.0, 104.0],
                              [-80.0,  -88.0, -96.0, -104.0],
                              [-180.0,  -188.0, -196.0, -204.0]], mask=False)
        SUN_AZ = np.ma.array([[148.0,  156.0, 164.0,  172.0],
                              [180.0,  188.0, 196.0, 204.0],
                              [180.0,  188.0, 196.0, 204.0],
                              [185.0,  193.0, 201.0, 209.0]], mask=False)

        RES = np.ma.array([[100., 100., 100., 100.],
                           [100., 100., 100., 100.],
                           [100.,  84.,  68.,  52.],
                           [5.,  21.,  37.,  53.]],
                          mask=False)
        rel_azi = get_absolute_azimuth_angle_diff(SAT_AZ, SUN_AZ)

        np.testing.assert_equal(rel_azi, RES)

def suite():
    """The suite for test_slerp
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestAngles))

    return mysuite


if __name__ == '__main__':
    unittest.main()
