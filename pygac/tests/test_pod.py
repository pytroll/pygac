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

from pygac.clock_offsets_converter import txt as clock_offsets_txt
from pygac.reader import ReaderError, NoTLEData
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
            1,  # 00...001
            QFlag.FATAL_FLAG,  # 100...00
            QFlag.CALIBRATION | QFlag.NO_EARTH_LOCATION,
            QFlag.TIME_ERROR | QFlag.DATA_GAP,
        ], dtype='>u4')
        # check if the bits look as expected
        bits = np.unpackbits(quality_indicators.view(np.uint8)).reshape((-1, 32))
        # For a big endian integer, the number 1 fills only the last of the 32 bits
        self.assertEqual(bits[0].sum(), 1)  # only one bit is filled
        self.assertEqual(bits[0][-1], 1)  # the last bit is filled
        # The fatal flag fills only the first bit
        self.assertEqual(bits[1].sum(), 1)  # only one bit is filled
        self.assertEqual(bits[1][0], 1)  # the first bit is filled

        # setup reader and test
        reader.scans = {self.reader._quality_indicators_key: quality_indicators}

        # default mask is QFlag.FATAL_FLAG | QFlag.CALIBRATION | QFlag.NO_EARTH_LOCATION
        expected_mask = np.array([False, True, True, False], dtype=bool)
        numpy.testing.assert_array_equal(reader.mask, expected_mask)

        # test individual flags
        expected_mask = np.array([False, False, False, True], dtype=bool)
        numpy.testing.assert_array_equal(
            reader._get_corrupt_mask(flags=QFlag.TIME_ERROR),
            expected_mask
        )
        # test combination of flags
        expected_mask = np.array([False, False, True, True], dtype=bool)
        flags = QFlag.DATA_GAP | QFlag.NO_EARTH_LOCATION
        numpy.testing.assert_array_equal(
            reader._get_corrupt_mask(flags=flags),
            expected_mask
        )

    @mock.patch('pygac.pod_reader.get_lonlatalt')
    @mock.patch('pygac.pod_reader.compute_pixels')
    @mock.patch('pygac.reader.Reader.get_tle_lines')
    @mock.patch('pygac.pod_reader.avhrr_gac')
    def test__adjust_clock_drift(self, avhrr_gac, get_tle_lines,
                                 compute_pixels, get_lonlatalt):
        """Test the clock drift adjustment."""
        sat_name = "fake_sat"
        reader = self.reader

        # We construct the following input
        # the scan lines do not start with 1 and have a gap
        scan_lines = np.array([15, 16, 17, 18, 22, 23, 24, 25])
        scan_rate = 0.5  # seconds/line
        # the first valid scans starts 1980-01-01 (without clock correction)
        # which leads to the following utcs for the given scan lines
        # ['1980-01-01T00:00:00.000', '1980-01-01T00:00:00.500',
        #  '1980-01-01T00:00:01.000', '1980-01-01T00:00:01.500',
        #  '1980-01-01T00:00:03.500', '1980-01-01T00:00:04.000',
        #  '1980-01-01T00:00:04.500', '1980-01-01T00:00:05.000']
        scan_utcs = (
            (1000 * scan_rate * (scan_lines - scan_lines[0])).astype('timedelta64[ms]')
            + np.datetime64("1980", "ms")
        )
        # For the geolocations, we assume an artificial swath of two pixels width
        # from north to south with constant lon, lat difference of 3deg for simplicity
        # lons = [[0, 3], [0, 3], [0, 3], [0, 3],
        #         [0, 3], [0, 3], [0, 3], [0, 3]],
        # lats = [[45, 45], [48, 48], [51, 51], [54, 54],
        #         [66, 66], [69, 69], [72, 72], [75, 75]]
        scan_angle = 3.0  # deg
        scan_lons, scan_lats = np.meshgrid(scan_angle*np.arange(2), scan_angle*scan_lines)

        # we assume a constant clock offset of 3.75 seconds
        # which should lead to the following adjustment on utcs
        # ['1979-12-31T23:59:56.250', '1979-12-31T23:59:56.750',
        #  '1979-12-31T23:59:57.250', '1979-12-31T23:59:57.750',
        #  '1979-12-31T23:59:59.750', '1980-01-01T00:00:00.250',
        #  '1980-01-01T00:00:00.750', '1980-01-01T00:00:01.250']
        offset = 3.75
        scan_offsets = offset*np.ones_like(scan_lines, dtype=float)  # seconds
        expected_utcs = scan_utcs - (1000*scan_offsets).astype('timedelta64[ms]')

        # the adjustment of geolocations should keep the lons unchanged,
        # but should shift the lats by scan_angel * offset / scan_rate
        # = 3deg/line * 3.75sec / 0.5sec/line = 22.5deg
        # [[22.5, 22.5], [25.5, 25.5], [28.5, 28.5], [31.5, 31.5],
        #  [43.5, 43.5], [46.5, 46.5], [49.5, 49.5], [52.5, 52.5]]
        expected_lons = scan_lons
        lats_shift = scan_angle*offset/scan_rate
        expected_lats = scan_lats - lats_shift

        # prepare the reader
        reader.scans = {"scan_line_number": scan_lines}
        reader.utcs = scan_utcs
        reader.lons = scan_lons
        reader.lats = scan_lats
        reader.spacecraft_name = sat_name

        # prepare offsets
        clock_offsets_txt[sat_name] = ("75001 000000 {offset}  "
                                       "85001 000000 {offset}").format(offset=offset)
        # set attitude coeffs
        reader._rpy = np.zeros(3)

        # set mocks for reader._compute_missing_lonlat call
        sgeom = mock.Mock()
        sgeom.times.return_value = [None]
        avhrr_gac.return_value = sgeom
        get_tle_lines.return_value = [None, None]
        compute_pixels.return_value = None
        # for the following mock, we need to know the missing values in advanced,
        # which we do, because we can convert the offset seconds into line number.
        offset_lines = offset / scan_rate
        min_line = np.floor(scan_lines[0] - offset_lines).astype(int)
        max_line = scan_lines[-1]
        missed_lines = np.setdiff1d(np.arange(min_line, max_line+1), scan_lines)
        n_missed = len(missed_lines)
        missed_lons = np.tile([0., scan_angle], n_missed)
        missed_lats = np.repeat(scan_angle*missed_lines, 2)
        get_lonlatalt.return_value = [missed_lons, missed_lats]

        # adjust clock drift
        reader._adjust_clock_drift()

        # check output
        # use allclose for geolocations, because the slerp interpolation
        # includes a transormation to cartesian coordinates and back to lon, lats.
        numpy.testing.assert_array_equal(reader.utcs, expected_utcs)
        numpy.testing.assert_allclose(reader.lons, expected_lons)
        numpy.testing.assert_allclose(reader.lats, expected_lats)

        # undo changes to clock_offsets_txt
        clock_offsets_txt.pop(sat_name)

    @mock.patch('pygac.pod_reader.get_offsets')
    @mock.patch('pygac.reader.Reader.get_tle_lines')
    def test__adjust_clock_drift_without_tle(self, get_tle_lines, get_offsets):
        """Test that clockdrift adjustment can handle missing TLE data."""
        reader = self.reader
        reader.utcs = np.zeros(10, dtype='datetime64[ms]')
        reader.scans = {"scan_line_number": np.arange(10)}
        get_offsets.return_value = np.zeros(10), np.zeros(10)
        get_tle_lines.side_effect = NoTLEData('No TLE data available')
        reader._adjust_clock_drift()  # should pass without errors
