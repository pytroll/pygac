#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author(s):

#   Stephan Finkensieper <stephan.finkensieper@dwd.de>

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

import unittest
import numpy as np
from pygac.gac_pod import GACPODReader


class TestPOD(unittest.TestCase):
    """Test the POD GAC reader"""

    longMessage = True

    def test_decode_timestamps(self):
        """Test POD timestamp decoding"""
        # Reference timestamps, one before 2000 one after 2000
        t2000_ref = (2001, 335, 53644260)
        t1900_ref = (1983, 336, 35058207)

        # Encoded version
        t2000_enc = np.array([847, 818, 35812])
        t1900_enc = np.array([42832, 534, 61983])

        # Test whether PODReader decodes them correctly
        self.assertEqual(GACPODReader.decode_timestamps(t2000_enc), t2000_ref,
                         msg='Timestamp after 2000 was decoded incorrectly')
        self.assertEqual(GACPODReader.decode_timestamps(t1900_enc), t1900_ref,
                         msg='Timestamp before 2000 was decoded incorrectly')


def suite():
    """The suite for test_pod"""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestPOD))
    return mysuite


if __name__ == '__main__':
    unittest.main()
