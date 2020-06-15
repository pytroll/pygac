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
"""Test the I/O."""

import unittest
try:
    import mock
except ImportError:
    from unittest import mock
import numpy as np
import numpy.testing
import pygac.gac_io as gac_io
import pygac.utils as utils


class TestIO(unittest.TestCase):
    """Test the gac_io module."""

    longMessage = True

    def test_strip_invalid_lat(self):
        """Test stripping the invalid lats."""
        lats = np.array([np.nan, 1, np.nan, 2, np.nan])
        start, end = utils.strip_invalid_lat(lats)
        self.assertEqual(start, 1)
        self.assertEqual(end, 3)

    def test_update_scanline(self):
        """Test updating the scanlines."""
        test_data = [{'new_start_line': 100, 'new_end_line': 200,
                      'scanline': 110, 'scanline_exp': 10},
                     {'new_start_line': 100, 'new_end_line': 200,
                      'scanline': 90, 'scanline_exp': None},
                     {'new_start_line': 100, 'new_end_line': 200,
                      'scanline': 210, 'scanline_exp': None}]
        for t in test_data:
            scanline_exp = t.pop('scanline_exp')
            scanline = utils._update_scanline(**t)
            self.assertEqual(scanline, scanline_exp)

    def test_slice(self):
        """Test slices."""
        ch = np.array([[1, 2, 3, 4, 5]]).transpose()
        sliced_exp = np.array([[2, 3, 4]]).transpose()

        # Without update
        sliced = utils._slice(ch, start_line=1, end_line=3)
        numpy.testing.assert_array_equal(sliced, sliced_exp)

        # With update
        sliced, updated = utils._slice(ch, start_line=1, end_line=3,
                                       update=[0, 2, 4, None])
        numpy.testing.assert_array_equal(sliced, sliced_exp)
        numpy.testing.assert_array_equal(updated, [None, 1, None, None])

        # Make sure slice is a copy
        ch += 1
        numpy.testing.assert_array_equal(sliced, sliced_exp)

    def test_slice_channel(self):
        """Test selection of user defined scanlines.

        Scanline Nr:  1  2  3  4  5  6  7  8  9 10 11 12 13 14
        Missing:      X        X                    X  X
        Inv. lat/lon:    X  X                             X  X

        Before stripping of invalid lats
        --------------------------------
        idx:             0  1     2  3  4  5  6  7        8  9
        data:            1  2     3  4  5  6  7  8        9 10

        After stripping of invalid lats
        -------------------------------
        idx:                      0  1  2  3  4  5
        data:                     3  4  5  6  7  8


        => All missing lines: [1, 2, 3, 4, 11, 12, 13, 14]
        """
        # Define test data
        ch = np.array([[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]).transpose()
        qualflags = np.array(
            [[2, 3, 5, 6, 7, 8, 9, 10, 13, 14]]).transpose()  # scanline number
        first_valid_lat = 2
        last_valid_lat = 7
        start_line = 1
        end_line = 4

        # Without lat stripping
        sliced_ref = np.array([[2, 3, 4, 5]]).transpose()
        sliced = utils.slice_channel(
            ch, start_line=start_line, end_line=end_line)
        numpy.testing.assert_array_equal(sliced, sliced_ref)

        # With lat stripping
        sliced_ref = np.array([[4, 5, 6, 7]]).transpose()
        sliced = utils.slice_channel(
            ch, start_line=start_line, end_line=end_line,
            first_valid_lat=first_valid_lat, last_valid_lat=last_valid_lat)
        numpy.testing.assert_array_equal(sliced, sliced_ref)

    def test_check_user_scanlines(self):
        """Check the scanlines."""
        # All scanlines
        start, end = utils.check_user_scanlines(0, 0, 100, 200)
        self.assertEqual(start, 0)
        self.assertEqual(end, 100)

        start, end = utils.check_user_scanlines(0, 0, None, None, 100)
        self.assertEqual(start, 0)
        self.assertEqual(end, 99)

        # Valid scanlines
        start, end = utils.check_user_scanlines(10, 20, 100, 200)
        self.assertEqual(start, 10)
        self.assertEqual(end, 20)

        start, end = utils.check_user_scanlines(10, 20, None, None, 100)
        self.assertEqual(start, 10)
        self.assertEqual(end, 20)

        # Invalid scanlines
        start, end = utils.check_user_scanlines(10, 110, 100, 200)
        self.assertEqual(start, 10)
        self.assertEqual(end, 100)

        start, end = utils.check_user_scanlines(10, 110, None, None, 100)
        self.assertEqual(start, 10)
        self.assertEqual(end, 99)

        self.assertRaises(ValueError, gac_io.check_user_scanlines,
                          110, 120, 100, 200)
        self.assertRaises(ValueError, gac_io.check_user_scanlines,
                          110, 120, None, None, 100)

    @mock.patch('pygac.gac_io.avhrrGAC_io')
    @mock.patch('pygac.gac_io.strip_invalid_lat')
    @mock.patch('pygac.gac_io.slice_channel')
    @mock.patch('pygac.gac_io.check_user_scanlines')
    def test_save_gac(self, check_user_scanlines, slice_channel,
                      strip_invalid_lat, avhrr_gac_io):
        """Test saving."""
        # Test scanline selection
        mm = mock.MagicMock()
        kwargs = dict(
            satellite_name=mm,
            xutcs=mm,
            lats=mm,
            lons=mm,
            ref1=mm,
            ref2=mm,
            ref3=mm,
            bt3=mm,
            bt4=mm,
            bt5=mm,
            sun_zen=mm,
            sat_zen=mm,
            sun_azi=mm,
            sat_azi=mm,
            rel_azi=mm,
            qual_flags=mm,
            gac_file=mm,
            meta_data=mm
        )
        slice_channel.return_value = mm
        strip_invalid_lat.return_value = 0, 0
        check_user_scanlines.return_value = 'start', 'end'

        gac_io.save_gac(start_line=0, end_line=0, **kwargs)
        slice_channel.assert_called_with(mock.ANY,
                                         start_line='start', end_line='end',
                                         first_valid_lat=mock.ANY,
                                         last_valid_lat=mock.ANY)
        expected_args = [
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
            mock.ANY,
            mock.ANY,
            mock.ANY,
            mock.ANY,
            mock.ANY,
            mock.ANY,
            mock.ANY,
            mock.ANY,
            'start',
            'end',
            mock.ANY,
            mock.ANY,
            mock.ANY,
            mock.ANY
        ]
        avhrr_gac_io.assert_called_with(*expected_args)


def suite():
    """Test suite for test_io."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestIO))
    return mysuite


if __name__ == '__main__':
    unittest.main()
