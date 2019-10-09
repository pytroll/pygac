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

import datetime
import unittest
import mock
import numpy as np
import numpy.testing
from pygac.gac_reader import GACReader


class TestGacReader(unittest.TestCase):
    """Test the common GAC Reader."""

    longMessage = True

    @mock.patch.multiple('pygac.gac_reader.GACReader', __abstractmethods__=set())
    def setUp(self, *mocks):
        self.reader = GACReader()

    def test_to_datetime64(self):
        """Test conversion from (year, jday, msec) to datetime64."""
        t0 = GACReader.to_datetime64(year=np.array(1970), jday=np.array(1),
                                     msec=np.array(0))
        self.assertEqual(t0.astype('i8'), 0,
                         msg='Conversion from (year, jday, msec) to datetime64 '
                             'is not correct')

    def test_to_datetime(self):
        """Test conversion from datetime64 to datetime."""
        dt = datetime.datetime(2019, 10, 1, 12, 15, 18, 12345)
        dt64 = np.datetime64(dt)
        self.assertEqual(self.reader.to_datetime(dt64), dt)

    def test_lineno2msec(self):
        """Test scanline timestamp estimation."""
        self.assertEqual(self.reader.lineno2msec(12345), 6172000)

    @mock.patch('pygac.gac_reader.GACReader._get_lonlat')
    @mock.patch('pygac.gac_reader.GACReader.get_corrupt_mask')
    @mock.patch('pygac.gac_reader.GACReader.adjust_clock_drift')
    @mock.patch('pygac.gac_reader.gtp.Gac_Lat_Lon_Interpolator')
    def test_get_lonlat(self, interpolator, adjust_clockdrift,
                        get_corrupt_mask, get_lonlat):
        """Test common lon/lat computation."""
        interpolator.return_value = np.array([np.nan, 1, 2, 3, -180.1, 180.1]), \
            np.array([1, 2, 3, np.nan, -90.1, 90.1])
        get_corrupt_mask.return_value = np.array([0, 0, 1, 0, 0, 0], dtype=bool)
        get_lonlat.return_value = None, None
        lons_exp = np.array([np.nan, 1, np.nan, 3., np.nan, np.nan])
        lats_exp = np.array([1, 2, np.nan, np.nan, np.nan, np.nan])
        lons, lats = self.reader.get_lonlat()
        get_lonlat.assert_called()
        adjust_clockdrift.assert_called()
        numpy.testing.assert_array_equal(lons, lons_exp)
        numpy.testing.assert_array_equal(lats, lats_exp)

        # Test caching
        methods = [get_lonlat, interpolator, adjust_clockdrift, get_corrupt_mask]
        for method in methods:
            method.reset_mock()
        self.reader.get_lonlat()
        for method in methods:
            method.asser_not_called()

    @mock.patch('pygac.gac_reader.GACReader._get_corrupt_mask')
    def test_get_corrupt_mask(self, get_corrupt_mask):
        """Test common computation of corrupt scanline mask."""
        get_corrupt_mask.return_value = [1, 2, 3]
        self.assertEqual(self.reader.get_corrupt_mask(), [1, 2, 3])

        # Test caching
        get_corrupt_mask.reset_mock()
        self.reader.get_corrupt_mask()
        get_corrupt_mask.assert_not_called()

    def test_midnight_scanline(self):
        """Test midnight scanline computation."""
        # Define test cases...
        # ... midnight scanline exists
        utcs1 = np.array([-3, -2, -1, 0, 1, 2, 3]).astype('datetime64[ms]')
        scanline1 = 2

        # ... midnight scanline does not exist
        utcs2 = np.array([1, 2, 3]).astype('datetime64[ms]')
        scanline2 = None

        for utcs, scanline in zip((utcs1, utcs2), (scanline1, scanline2)):
            self.reader.utcs = utcs
            self.reader.times = utcs.astype(datetime.datetime)
            self.assertEqual(self.reader.get_midnight_scanline(), scanline,
                             msg='Incorrect midnight scanline')

    def test_miss_lines(self):
        """Test detection of missing scanlines."""
        lines = [2, 4, 5, 6, 10, 11, 12]
        miss_lines_ref = [1, 3, 7, 8, 9]
        self.reader.scans = np.zeros(len(lines), dtype=[('scan_line_number', 'i2')])
        self.reader.scans['scan_line_number'] = lines
        self.assertTrue((self.reader.get_miss_lines() == miss_lines_ref).all(),
                        msg='Missing scanlines not detected correctly')

    def test_tle2datetime64(self, *mocks):
        """Test conversion from TLE timestamps to datetime64."""
        dates = np.array([70365.1234, 18001.25])
        dates64_exp = [np.datetime64(datetime.datetime(1970, 12, 31, 2, 57, 41, 760000), '[ms]'),
                       np.datetime64(datetime.datetime(2018, 1, 1, 6, 0), '[ms]')]
        dates64 = GACReader.tle2datetime64(dates)
        self.assertTrue(np.all(dates64 == dates64_exp))

    @mock.patch('pygac.gac_reader.GACReader.get_times')
    @mock.patch('pygac.gac_reader.GACReader.get_tle_file')
    def test_get_tle_lines(self, get_tle_file, *mocks):
        """Test identification of closest TLE lines."""
        tle_data = ['1 38771U 12049A   18363.63219793 -.00000013  00000-0  14176-4 0  9991\r\n',
                    '2 38771  98.7297  60.1350 0002062  95.9284  25.0713 14.21477560325906\r\n',
                    '1 38771U 12049A   18364.62426010 -.00000015  00000-0  13136-4 0  9990\r\n',  # 2018-12-30 14:58
                    '2 38771  98.7295  61.1159 0002062  94.5796  60.2561 14.21477636326047\r\n',
                    '1 38771U 12049A   18364.94649306 -.00000018  00000-0  12040-4 0  9996\r\n',  # 2018-12-30 22:42
                    '2 38771  98.7295  61.4345 0002060  94.1226 268.7521 14.21477633326092\r\n',
                    '1 38771U 12049A   18365.81382142 -.00000015  00000-0  13273-4 0  9991\r\n',
                    '2 38771  98.7294  62.2921 0002057  92.7653  26.0030 14.21477711326215\r\n']

        expected = {
            datetime.datetime(2018, 12, 20, 12, 0): None,
            datetime.datetime(2018, 12, 28, 12, 0): 0,
            datetime.datetime(2018, 12, 30, 16, 0): 2,
            datetime.datetime(2018, 12, 30, 20, 0): 4,
            datetime.datetime(2019, 1, 1, 12, 0): 6,
            datetime.datetime(2019, 1, 8, 12, 0): None
        }

        get_tle_file.return_value = tle_data
        for time, tle_idx in expected.items():
            self.reader.times = [time]
            self.reader.tle_lines = None
            if tle_idx is None:
                self.assertRaises(IndexError, self.reader.get_tle_lines)
            else:
                tle1, tle2 = self.reader.get_tle_lines()
                self.assertEqual(tle1, tle_data[tle_idx])
                self.assertEqual(tle2, tle_data[tle_idx + 1])

    @mock.patch('pygac.gac_reader.GACReader.get_tsm_pixels')
    def test_mask_tsm_pixels(self, get_tsm_pixels):
        """Test masking of pixels affected by the scan motor issue."""
        get_tsm_pixels.return_value = ([0, 1], [0, 1])
        channels = np.array([[[1., 2., 3.],
                              [1., 2., 3.]],
                             [[1., 2., 3.],
                              [1., 2., 3.]]])  # (lines, pixels, channels)
        masked_exp = np.array([[[np.nan, np.nan, np.nan],
                                [1., 2., 3.]],
                               [[1., 2., 3.],
                                [np.nan, np.nan, np.nan]]])
        self.reader.mask_tsm_pixels(channels)  # masks in-place
        numpy.testing.assert_array_equal(channels, masked_exp)

    def _get_scanline_numbers(self):
        """Create artificial scanline numbers with some corruptions.

        Returns:
            Corrupted and corrected scanline numbers.
        """
        along_track = 12000
        scans = np.zeros(12000, dtype=[("scan_line_number", ">u2")])
        scans["scan_line_number"] = np.arange(1, along_track+1)

        # ... with 500 missing scanlines at scanline 8000
        scans["scan_line_number"][8000:] += 500
        corrected = scans["scan_line_number"].copy()

        # ... and some spikes here and there
        scans["scan_line_number"][3000] += 1E4
        scans["scan_line_number"][9000] -= 1E4
        corrected = np.delete(corrected, [3000, 9000])

        return scans, corrected

    def test_correct_scan_line_numbers(self):
        """Test scanline number correction."""
        scans, expected = self._get_scanline_numbers()
        self.reader.scans = scans
        self.reader.correct_scan_line_numbers()
        numpy.testing.assert_array_equal(self.reader.scans['scan_line_number'],
                                         expected)

    @mock.patch('pygac.gac_reader.GACReader.get_header_timestamp')
    def test_correct_times_thresh(self, get_header_timestamp):
        """Test correction of scanline timestamps."""
        header_time = datetime.datetime(2016, 8, 16, 16, 7, 36)

        # Create artificial timestamps
        _, scan_line_numbers = self._get_scanline_numbers()
        t0 = np.array([header_time], dtype="datetime64[ms]").astype("i8")[0]
        shift = 1000
        msecs = t0 + shift + scan_line_numbers / GACReader.scan_freq
        utcs_expected = msecs.copy().astype(">M8[ms]")

        # Add some corruptions
        msecs[3000:] += 1800 * 1000
        msecs[1000] += 24*3600*1000
        msecs[2000] -= 24*3600*1000

        # Mock reader
        get_header_timestamp.return_value = header_time
        self.reader.utcs = msecs.astype(">M8[ms]")
        self.reader.scans = np.array(scan_line_numbers,
                                     dtype=[("scan_line_number", ">u2")])

        # Test correction
        self.reader.correct_times_thresh()
        numpy.testing.assert_array_equal(self.reader.utcs, utcs_expected)


def suite():
    """The suite for test_reader"""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestGacReader))
    return mysuite


if __name__ == '__main__':
    unittest.main()
