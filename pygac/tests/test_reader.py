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
from pygac.gac_reader import GACReader


class TestGacReader(unittest.TestCase):
    """Test the common GAC Reader"""

    longMessage = True

    def test_to_datetime64(self):
        """Test conversion from (year, jday, msec) to datetime64"""
        t0 = GACReader.to_datetime64(year=np.array(1970), jday=np.array(1),
                                     msec=np.array(0))
        self.assertEqual(t0.astype('i8'), 0,
                         msg='Conversion from (year, jday, msec) to datetime64 '
                             'is not correct')

    @mock.patch.multiple('pygac.gac_reader.GACReader', __abstractmethods__=set())
    def test_midnight_scanline(self, *mocks):
        """Test midnight scanline computation"""
        reader = GACReader()

        # Define test cases...
        # ... midnight scanline exists
        utcs1 = np.array([-3, -2, -1, 0, 1, 2, 3]).astype('datetime64[ms]')
        scanline1 = 2

        # ... midnight scanline does not exist
        utcs2 = np.array([1, 2, 3]).astype('datetime64[ms]')
        scanline2 = None

        for utcs, scanline in zip((utcs1, utcs2), (scanline1, scanline2)):
            reader.utcs = utcs
            reader.times = utcs.astype(datetime.datetime)
            self.assertEqual(reader.get_midnight_scanline(), scanline,
                             msg='Incorrect midnight scanline')

    @mock.patch.multiple('pygac.gac_reader.GACReader', __abstractmethods__=set())
    def test_miss_lines(self, *mocks):
        """Test detection of missing scanlines"""
        reader = GACReader()
        lines = [2, 4, 5, 6, 10, 11, 12]
        miss_lines_ref = [1, 3, 7, 8, 9]
        reader.scans = np.zeros(len(lines), dtype=[('scan_line_number', 'i2')])
        reader.scans['scan_line_number'] = lines
        self.assertTrue((reader.get_miss_lines() == miss_lines_ref).all(),
                        msg='Missing scanlines not detected correctly')

    def test_tle2datetime64(self, *mocks):
        """Test conversion from TLE timestamps to datetime64"""
        dates = np.array([70365.1234, 18001.25])
        dates64_exp = [np.datetime64(datetime.datetime(1970, 12, 31, 2, 57, 41, 760000), '[ms]'),
                       np.datetime64(datetime.datetime(2018, 1, 1, 6, 0), '[ms]')]
        dates64 = GACReader.tle2datetime64(dates)
        self.assertTrue(np.all(dates64 == dates64_exp))

    @mock.patch.multiple('pygac.gac_reader.GACReader', __abstractmethods__=set())
    @mock.patch('pygac.gac_reader.GACReader.get_times')
    @mock.patch('pygac.gac_reader.GACReader.get_tle_file')
    def test_get_tle_lines(self, get_tle_file, *mocks):
        """Test identification of closest TLE lines"""
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
        reader = GACReader(tle_thresh=7)
        for time, tle_idx in expected.items():
            reader.times = [time]
            reader.tle_lines = None
            if tle_idx is None:
                self.assertRaises(IndexError, reader.get_tle_lines)
            else:
                tle1, tle2 = reader.get_tle_lines()
                self.assertEqual(tle1, tle_data[tle_idx])
                self.assertEqual(tle2, tle_data[tle_idx + 1])


def suite():
    """The suite for test_reader"""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestGacReader))
    return mysuite


if __name__ == '__main__':
    unittest.main()
