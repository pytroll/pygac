#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author(s):

#   Stephan Finkensieper <stephan.finkensieper@dwd.de>
#   Carlos Horn <carlos.horn@external.eumetsat.int>

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
"""Test module for the pod reading."""

import datetime as dt
import unittest
import numpy as np
import numpy.testing
import sys
try:
    from unittest import mock
except ImportError:
    import mock

from pygac.reader import ReaderError
from pygac.gac_pod import GACPODReader
from pygac.lac_pod import LACPODReader
from pygac.tests.utils import CalledWithArray


class TestPOD(unittest.TestCase):
    """Test the POD GAC reader."""

    longMessage = True

    def setUp(self):
        """Set up the test."""
        self.reader = GACPODReader()
        # python 2 compatibility
        if sys.version_info.major < 3:
            self.assertRaisesRegex = self.assertRaisesRegexp

    def test__validate_header(self):
        """Test the header validation"""
        filename = b'NSS.GHRR.TN.D80001.S0332.E0526.B0627173.WI'
        head = {'data_set_name': filename}
        GACPODReader._validate_header(head)
        # wrong name pattern
        with self.assertRaisesRegex(ReaderError,
                                    'Data set name .* does not match!'):
            head = {'data_set_name': b'abc.txt'}
            GACPODReader._validate_header(head)
        # wrong platform
        name = b'NSS.GHRR.NL.D02187.S1904.E2058.B0921517.GC'
        with self.assertRaisesRegex(ReaderError,
                                    'Improper platform id "NL"!'):
            head = {'data_set_name': name}
            GACPODReader._validate_header(head)
        # wrong transfer mode
        name = filename.replace(b'GHRR', b'LHRR')
        with self.assertRaisesRegex(ReaderError,
                                    'Improper transfer mode "LHRR"!'):
            head = {'data_set_name': name}
            GACPODReader._validate_header(head)
        # other change reader
        head = {'data_set_name': name}
        LACPODReader._validate_header(head)

    @mock.patch('pygac.reader.Reader.get_calibrated_channels')
    def test__get_calibrated_channels_uniform_shape(self, get_channels):
        """Test the uniform shape as required by gac_io.save_gac."""
        channels = np.arange(2*2*5, dtype=float).reshape((2, 2, 5))
        get_channels.return_value = channels
        uniform_channels = self.reader._get_calibrated_channels_uniform_shape()
        self.assertTrue(np.isnan(uniform_channels[:, :, 2]).all())
        self.assertTrue(uniform_channels[:, :, [0, 1, 3, 4, 5]].sum()
                        == channels.sum())

    def test_decode_timestamps(self):
        """Test POD timestamp decoding."""
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
        """Test getting times."""
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

    @mock.patch('pygac.pod_reader.get_tsm_idx')
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

    def test_quality_indicators(self):
        """Test the quality indicator unpacking."""
        reader = self.reader
        QFlag = reader.QFlag
        quality_indicators = np.array([
            0,  # nothing flagged
            -1,  # everything flagged
            QFlag.CALIBRATION | QFlag.NO_EARTH_LOCATION,
            QFlag.TIME_ERROR | QFlag.DATA_GAP,
            QFlag.FATAL_FLAG
        ], dtype=np.uint32)
        reader.scans = {self.reader._quality_indicators_key: quality_indicators}
        # test mask, i.e. QFlag.FATAL_FLAG | QFlag.CALIBRATION | QFlag.NO_EARTH_LOCATION
        expected_mask = np.array([False, True, True, False, True], dtype=bool)
        numpy.testing.assert_array_equal(reader.mask, expected_mask)
        # test individual flags
        self.assertTrue(reader._get_corrupt_mask(flags=QFlag.FATAL_FLAG).any())
        # count the occurence (everything flagged and last entrance => 2)
        self.assertEqual(reader._get_corrupt_mask(flags=QFlag.FATAL_FLAG).sum(), 2)


def suite():
    """Test suite for test_pod."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestPOD))
    return mysuite


if __name__ == '__main__':
    unittest.main()
