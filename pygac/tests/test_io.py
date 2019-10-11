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
import mock
import numpy as np
import numpy.testing
import pygac.gac_io as gac_io
import pygac.utils as utils
from pygac.tests.utils import CalledWithArray


class TestIO(unittest.TestCase):
    """Test the gac_io module"""

    longMessage = True

    def test_strip_invalid_lat(self):
        lats = np.array([np.nan, 1, np.nan, 2, np.nan])
        start, end = utils.strip_invalid_lat(lats)
        self.assertEqual(start, 1)
        self.assertEqual(end, 3)

    def test_update_start_end_line(self):
        test_data = [{'user_start': 100, 'user_end': 200,
                      'valid_lat_start': 0, 'valid_lat_end': 300,
                      'start_exp': 100, 'end_exp': 200},
                     {'user_start': 100, 'user_end': 200,
                      'valid_lat_start': 0, 'valid_lat_end': 175,
                      'start_exp': 100, 'end_exp': 175},
                     {'user_start': 100, 'user_end': 200,
                      'valid_lat_start': 125, 'valid_lat_end': 300,
                      'start_exp': 125, 'end_exp': 200},
                     {'user_start': 100, 'user_end': 200,
                      'valid_lat_start': 125, 'valid_lat_end': 175,
                      'start_exp': 125, 'end_exp': 175}]
        for t in test_data:
            start_exp = t.pop('start_exp')
            end_exp = t.pop('end_exp')
            start, end = utils.update_start_end_line(**t)
            self.assertEqual(start, start_exp)
            self.assertEqual(end, end_exp)

    def test_update_scanline(self):
        test_data = [{'new_start_line': 100, 'new_end_line': 200,
                      'scanline': 110, 'scanline_exp': 10},
                     {'new_start_line': 100, 'new_end_line': 200,
                      'scanline': 90, 'scanline_exp': None},
                     {'new_start_line': 100, 'new_end_line': 200,
                      'scanline': 210, 'scanline_exp': None}]
        for t in test_data:
            scanline_exp = t.pop('scanline_exp')
            scanline = utils.update_scanline(**t)
            self.assertEqual(scanline, scanline_exp)

    def test_update_missing_scanlines(self):
        qual_flags = np.array([[1, 2, 4, 5, 6, 8, 9, 11, 12]]).transpose()
        miss_lines = np.array([3, 7, 10])
        test_data = [{'valid_lat_start': 0, 'valid_lat_end': 8,
                      'miss_lines_exp': [3, 7, 10]},
                     {'valid_lat_start': 3, 'valid_lat_end': 6,
                      'miss_lines_exp': [1, 2, 3, 4, 7, 10, 11, 12]}]
        for t in test_data:
            miss_lines_exp = t.pop('miss_lines_exp')
            miss_lines = utils.update_missing_scanlines(
                miss_lines=miss_lines, qual_flags=qual_flags, **t)
            numpy.testing.assert_array_equal(miss_lines, miss_lines_exp)

    @mock.patch('pygac.utils.update_scanline')
    @mock.patch('pygac.utils.update_missing_scanlines')
    @mock.patch('pygac.utils.update_start_end_line')
    def test_slice(self, update_start_end_line, update_missing_scanlines,
                   update_scanline):
        ch = np.array([[1, 2, 3, 4, 5]]).transpose()
        start_line = 'start_line'
        end_line = 'end_line'
        valid_lat_start = 1
        valid_lat_end = 3

        update_start_end_line.return_value = 1, 3
        update_missing_scanlines.return_value = 'miss_lines'
        update_scanline.return_value = 'midnight'

        sliced_exp = np.array([[2, 3, 4]]).transpose()

        # Test slicing without lat stripping and other scanline updates
        sliced = utils.slice_channel(ch,
                                     start_line=start_line,
                                     end_line=end_line)
        numpy.testing.assert_array_equal(sliced, sliced_exp)
        update_start_end_line.assert_called_with(
            user_start=start_line,
            user_end=end_line,
            valid_lat_start=0,
            valid_lat_end=4)
        update_missing_scanlines.assert_not_called()
        update_scanline.assert_not_called()

        # Test slicing with lat stripping and other scanline updates
        sliced, miss_lines, midnight_scanline = utils.slice_channel(
            ch,
            start_line=start_line,
            end_line=end_line,
            valid_lat_start=valid_lat_start,
            valid_lat_end=valid_lat_end,
            miss_lines='dummy',
            midnight_scanline='dummy',
            qual_flags='dummy')
        numpy.testing.assert_array_equal(sliced, sliced_exp)
        update_start_end_line.assert_called_with(
            user_start=start_line,
            user_end=end_line,
            valid_lat_start=valid_lat_start,
            valid_lat_end=valid_lat_end)
        update_missing_scanlines.assert_called_with(
            miss_lines='dummy',
            qual_flags='dummy',
            valid_lat_start=valid_lat_start,
            valid_lat_end=valid_lat_end)
        update_scanline.assert_called_with('dummy',
                                           new_start_line=1,
                                           new_end_line=3)

        # Make sure slice is a copy
        ch += 1
        numpy.testing.assert_array_equal(sliced, sliced_exp)

    def test_save_gac(self):
        """Test selection of user defined scanlines

        Scanline Nr:  1  2  3  4  5  6  7  8  9 10 11 12 13 14
        Missing:      X           X                 X        X
        Inv. lat/lon:    X  X  X                 X     X  X

        => All missing lines: [1, 2, 3, 4, 5, 10, 11, 12, 13, 14
        """
        # Define test data
        along = 10
        across = 3
        miss_lines = np.array([1, 5, 11, 14])  # never recorded
        qualflags = np.array(
            [[2, 3, 4, 6, 7, 8, 9, 10, 12, 13]]).transpose()  # scanline number
        lats = np.arange(along*across, dtype=float).reshape((along, across))
        lats[:2+1, :] = np.nan  # Make first 3 lines invalid
        lats[7:, :] = np.nan  # Make last 3 lines invalid
        lons = lats.copy()
        xutcs = np.arange(along).astype('datetime64[ms]')
        midnight_scanline = 5
        dummydata = np.ones(lats.shape)
        start_end_lines = [(0, 0), (4, 5), (0, 4), (6, 9)]

        # Define reference data for each pair of start/end lines. Lines with
        # invalid lat/lon info (here 0, 1, 2, 7, 8, 9) must be ignored
        new_start_lines = [3, 4, 3, 6]
        new_end_lines = [6, 5, 4, 6]
        midnight_scanlines_ref = [2, 1, None, None]
        all_miss_lines_ref = np.array([1, 2, 3, 4, 5, 10, 11, 12, 13, 14])

        for itest, (start_line, end_line) in enumerate(start_end_lines):
            new_start_line = new_start_lines[itest]
            new_end_line = new_end_lines[itest]
            lats_ref = 1000.0 * lats[new_start_line:new_end_line + 1, :]
            lats_ref[np.isnan(lats_ref)] = gac_io.MISSING_DATA_LATLON
            qualflags_ref = qualflags[new_start_line:new_end_line+1]
            xutcs_ref = xutcs[new_start_line:new_end_line + 1]
            midnight_scanline_ref = midnight_scanlines_ref[itest]

            with mock.patch('pygac.gac_io.avhrrGAC_io') as io_mock:
                gac_io.save_gac(
                    satellite_name='dummy',
                    xutcs=xutcs,
                    lats=lats.copy(),
                    lons=lons.copy(),
                    ref1=dummydata,
                    ref2=dummydata,
                    ref3=dummydata,
                    bt3=dummydata,
                    bt4=dummydata,
                    bt5=dummydata,
                    sun_zen=dummydata,
                    sat_zen=dummydata,
                    sun_azi=dummydata,
                    sat_azi=dummydata,
                    rel_azi=dummydata,
                    qual_flags=qualflags,
                    start_line=start_line,
                    end_line=end_line,
                    gac_file='dummy',
                    midnight_scanline=midnight_scanline,
                    miss_lines=miss_lines,
                )
                expected_args = [
                    mock.ANY,
                    CalledWithArray(xutcs_ref),
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    CalledWithArray(lats_ref),
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    CalledWithArray(qualflags_ref),
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    midnight_scanline_ref,
                    CalledWithArray(all_miss_lines_ref)
                ]
                io_mock.assert_called_with(*expected_args)


def suite():
    """The suite for test_io"""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestIO))
    return mysuite


if __name__ == '__main__':
    unittest.main()
