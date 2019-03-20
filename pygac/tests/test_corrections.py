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

"""Test corrections of scanline numbers, timestamps and scan motor issue"""


import unittest
from pygac.gac_reader import GACReader
import pygac.correct_tsm_issue as tsm

import numpy as np
import datetime


TEST_TIMESTAMP = datetime.datetime(2016, 8, 16, 16, 7, 36)
TEST_SPACECRAFT_ID = 1


class TestGACReader(GACReader):
    tsm_affected_intervals = {
        TEST_SPACECRAFT_ID: [(TEST_TIMESTAMP.replace(hour=0, minute=0),
                             TEST_TIMESTAMP.replace(hour=23, minute=59))]
    }

    def read(self, filename):
        pass

    def get_header_timestamp(self):
        return TEST_TIMESTAMP

    def _get_times(self):
        pass


class TestCorrections(unittest.TestCase):
    """Test corrections of scanline numbers, timestamps and scan motor issue"""

    def setUp(self):
        # Create artificial scanline numbers
        along_track = 12000
        scans = np.zeros(12000, dtype=[("scan_line_number", ">u2")])
        scans["scan_line_number"] = np.arange(1, along_track+1)

        # ... with 200 missing scanlines at scanline 8000
        scans["scan_line_number"][8000:] += 500
        self.correct_scan_line_numbers = scans["scan_line_number"].copy()

        # ... and some spikes here and there
        scans["scan_line_number"][3000] += 1E4
        scans["scan_line_number"][9000] -= 1E4
        self.correct_scan_line_numbers = np.delete(
            self.correct_scan_line_numbers, [3000, 9000])

        # Create artificial timestamps
        scanning_rate = 2.0 / 1000.0  # scanlines per millisecond
        t0 = np.array([TEST_TIMESTAMP.isoformat()],
                      dtype="datetime64[ms]").astype("i8")[0]
        shift = 1000
        msecs = t0 + shift + self.correct_scan_line_numbers/scanning_rate
        self.correct_utcs = msecs.copy().astype(">M8[ms]")

        # ... with some corruptions here and there
        msecs[1000] += 24*3600*1000
        msecs[2000] -= 24*3600*1000
        msecs[3000:] += 1800*1000

        # Create a GAC Reader and parse the test data
        self.reader = TestGACReader()
        self.reader.filename = 'test'
        self.reader.scans = scans
        self.reader.utcs = msecs.astype(">M8[ms]")
        self.reader.times = self.reader.utcs.astype(datetime.datetime)

        # Create artificial channel data with some noisy regions. Channel order:
        # 1, 2, 3b, 4, 5, 3a
        cols = np.arange(409)
        rows = np.arange(12000)
        mcols, mrows = np.meshgrid(cols, rows)
        data = np.sin(mrows/float(rows.size))*np.cos(mcols/float(cols.size))
        channels = [data.copy() for i in range(6)]

        # ... add noise to channel 1 and 4
        channels_noise = [c.copy() for c in channels]
        noise_cols, noise_rows = np.meshgrid(np.arange(200, 300),
                                             np.arange(4000, 5000))
        ampl = 0.5
        for ichan in (0, 3):
            channels_noise[ichan][noise_rows, noise_cols] += ampl*np.random.rand(
                noise_rows.shape[0], noise_rows.shape[1])

        # ... scale data
        for i in (0, 1, 5):
            # Reflectances
            channels[i] /= 0.01*0.01
            channels_noise[i] /= 0.01 * 0.01
        for i in (2, 3, 4):
            # Brightness temperatures
            channels[i] = (channels[i] - 273.15) / 0.01
            channels_noise[i] = (channels_noise[i] - 273.15)/0.01

        self.channels = channels
        self.channels_noise = channels_noise
        self.noise_cols = noise_cols
        self.noise_rows = noise_rows

    def test_scan_line_number(self):
        """Test scan line number correction"""
        self.reader.correct_scan_line_numbers()
        self.assertTrue(np.allclose(self.correct_scan_line_numbers,
                                    self.reader.scans["scan_line_number"]),
                        msg='Scanline number correction failed')

    def test_timestamp(self):
        """Test timestamp correction"""
        self.reader.correct_scan_line_numbers()
        self.reader.correct_times_thresh()
        self.assertTrue(np.allclose(self.correct_utcs.astype('i8'),
                                    self.reader.utcs.astype('i8')),
                        msg='Timestamp correction failed')

    def test_scan_motor_issue(self):
        """Test correction of the scan motor issue"""
        fv = -32001
        ch1, ch2, ch3b, ch4, ch5, ch3a = self.channels_noise
        channels_corr = tsm.flag_pixels(channel1=ch1, channel2=ch2,
                                        channel3b=ch3b, channel4=ch4,
                                        channel5=ch5, channel3a=ch3a,
                                        fillv=fv)
        for channel, channel_corr in zip(self.channels, channels_corr):
            nonoise = channel_corr != fv
            self.assertTrue(
                np.allclose(channel_corr[self.noise_rows, self.noise_cols], fv),
                msg='Some corrupt pixels have not been masked')
            self.assertTrue(
                np.allclose(channel_corr[nonoise] - channel[nonoise], 0.0),
                msg='Non-corrupt pixels have been modified')

    def test_is_tsm_affected(self):
        """Test identification of TSM affected orbits"""
        # Affected platform and time interval
        self.reader.spacecraft_id = TEST_SPACECRAFT_ID
        self.assertTrue(self.reader.is_tsm_affected(),
                        msg='Affected orbit is not recognized')

        # Unaffected platform
        self.reader.spacecraft_id = TEST_SPACECRAFT_ID + 31415
        self.assertFalse(self.reader.is_tsm_affected(),
                         msg='Unaffected platform is recognized mistakenly')

        # Unaffected time
        self.reader.spacecraft_id = TEST_SPACECRAFT_ID
        self.reader.times = [t + datetime.timedelta(days=365*12)
                             for t in self.reader.times]
        self.assertFalse(self.reader.is_tsm_affected(),
                         msg='Unaffected time is recognized mistakenly')


class MeanFilterTest(unittest.TestCase):
    """Test gridbox mean filter"""

    def runTest(self):
        # Define test data
        data = np.ma.array(
            [[1, 2, 2, 1],
             [2, 1, 2, 1],
             [1, 1, 2, 2],
             [2, 2, 1, 1]]
        )
        data.mask = np.zeros(data.shape)
        data.mask[1, 1:3] = 1
        data.mask[3, -1] = 1

        # Define reference
        filtered_ref = np.array(
            [[5/3., 7/4., 6/4., 4/3.],
             [7/5., 11/7., 11/7., 8/5.],
             [8/5., 11/7., 9/6., 6/4.],
             [6/4., 9/6., 8/5., 5/3.]]
        )

        # Apply mean
        filtered = tsm.mean_filter(data=data, box_size=3, fill_value=-999)

        # Compare results against reference
        self.assertTrue(np.allclose(filtered, filtered_ref),
                        msg='Mean filter produces incorrect results.')


def suite():
    """The suite for test_corrections"""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestCorrections))
    mysuite.addTest(loader.loadTestsFromTestCase(MeanFilterTest))

    return mysuite


if __name__ == '__main__':
    unittest.main()
