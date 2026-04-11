#!/usr/bin/env python

# Copyright (c) 2014-2025 Pytroll Developers

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

"""Test function for the calibration PRT numbering for NOAA calibration
"""

import sys
import unittest

import numpy as np

from pygac.calibration.noaa import get_prt_nos

#
# dummy user data of line number and prt values (average of three values per
# scanline) to test allocation of which PRT number we are at. Includes
# glitches.
# Data taken from a section of real data 
class TestPrtNumbering(unittest.TestCase):

    def test_prt_numbering(self):
        #
        # Data taken from real problematic case
        #
        line_number = np.array([7386,7387,7388,7389,7390,7391,
                                7392,7393,7394,7395,
                                7396,7397,7400,7401,7402,
                                7403,7404,7405,7406,7407,
                                7408,7409])
        prt_counts = np.array([0.0,270.0,263.0,269.33,270.0,
                               0.0,269.0,263.0,269.33,
                               270.0,0.0,0.0,0.0,263.0,
                               269.0,270.0,0.0,270.0,263.0,
                               269.0,270.0,0.0])
        prt_threshold = 50
        gac = True
        
        expected_iprt = [0,3,1,4,2,0,3,1,4,2,0,0,0,3,1,4,0,3,1,4,2,0]

        iprt = get_prt_nos(prt_counts,prt_threshold,line_number,gac)

        np.testing.assert_allclose(expected_iprt,iprt)
        
