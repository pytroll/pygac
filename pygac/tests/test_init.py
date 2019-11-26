#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014-2019 Pytroll Developers

# Author(s):

#   Nina Hakansson <nina.hakansson@smhi.se>

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

"""Test function for the angle calculation."""

import unittest

import numpy as np

from pygac import calculate_sun_earth_distance_correction


class TestInit(unittest.TestCase):
    """Test function for the angle calculation."""

    def test_calculate_sun_earth_distance_correction(self):
        """Test function for the sun distance corretction."""
        corr = calculate_sun_earth_distance_correction(3)
        np.testing.assert_almost_equal(corr, 0.96660494, decimal=7)


def suite():
    """Test non angle functions in pygac init file."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestInit))

    return mysuite


if __name__ == '__main__':
    unittest.main()
