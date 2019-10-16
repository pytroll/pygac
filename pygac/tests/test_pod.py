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

import datetime as dt
import unittest
import numpy as np
import numpy.testing
try:
    from unittest import mock
except ImportError:
    import mock

from pygac.gac_pod import GACPODReader
from pygac.tests.utils import CalledWithArray


class TestPOD(unittest.TestCase):
    """Test the POD GAC reader"""

    longMessage = True

    def setUp(self):
        self.reader = GACPODReader()

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

    @mock.patch('pygac.gac_pod.GACPODReader.decode_timestamps')
    def test_get_header_timestamp(self, decode_timestamps):
        """Test readout of header timestamp."""
        self.reader.head = {'start_time': 123}
        decode_timestamps.return_value = np.array(
            [2019]), np.array([123]), np.array([123456])
        time = self.reader.get_header_timestamp()
        decode_timestamps.assert_called_with(123)
        self.assertEqual(time, dt.datetime(2019, 5, 3, 0, 2, 3, 456000))

    @mock.patch('pygac.gac_pod.GACPODReader.decode_timestamps')
    def test_get_times(self, decode_timestamps):
        self.reader.scans = {'time_code': 123}
        self.reader._get_times()
        decode_timestamps.assert_called_with(123)

    def test_get_lonlat(self):
        """Test readout of lon/lat coordinates."""
        earth_loc = 128 * np.array([[1, 2, 3, 4],
                                    [5, 6, 7, 8]])
        self.reader.scans = {'earth_location': earth_loc}

        lons_exp = np.array([[2, 4],
                             [6, 8]])
        lats_exp = np.array([[1, 3],
                             [5, 7]])

        lons, lats = self.reader._get_lonlat()
        numpy.testing.assert_array_equal(lons, lons_exp)
        numpy.testing.assert_array_equal(lats, lats_exp)

    @mock.patch('pygac.gac_pod.get_tsm_idx')
    def test_get_tsm_pixels(self, get_tsm_idx):
        """Test channel set used for TSM correction."""
        ones = np.ones((409, 100))
        zeros = np.zeros(ones.shape)
        ch1 = 1*ones
        ch2 = 2*ones
        ch4 = 4*ones
        ch5 = 5*ones
        channels = np.dstack((ch1, ch2, zeros, ch4, ch5))
        self.reader.get_tsm_pixels(channels)
        get_tsm_idx.assert_called_with(CalledWithArray(ch1),
                                       CalledWithArray(ch2),
                                       CalledWithArray(ch4),
                                       CalledWithArray(ch5))


def suite():
    """The suite for test_pod"""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestPOD))
    return mysuite


if __name__ == '__main__':
    unittest.main()
