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
from pygac.gac_io import update_start_end_line


class CalledWithAdapter(object):
    def __init__(self, array):
        self.array = array

    def __repr__(self):
        return repr(self.array)

    def __eq__(self, other):
        return np.all(self.array == other)


class TestIO(unittest.TestCase):
    """Test the gac_io module"""

    longMessage = True

    def test_update_start_end(self):
        x = np.arange(0, 20)
        start_line = 5
        end_line = 15

        for temp_start_line in (0, start_line, start_line + 1):
            for temp_end_line in (x.size, end_line, end_line - 1):
                where = 'temporary start/end lines {0}/{1}'.format(
                    temp_start_line, temp_end_line)

                # Define reference
                if temp_start_line > start_line:
                    ref_start = temp_start_line
                else:
                    ref_start = start_line
                if temp_end_line < end_line:
                    ref_end = temp_end_line
                else:
                    ref_end = end_line
                ref = np.arange(ref_start, ref_end + 1)

                # Compute new start & end line
                new_start_line, new_end_line = update_start_end_line(
                    start_line=start_line,
                    end_line=end_line,
                    temp_start_line=temp_start_line,
                    temp_end_line=temp_end_line)

                # Slice twice (as in gac_io)
                slice1 = x[temp_start_line:temp_end_line + 1].copy()
                self.assertGreaterEqual(
                    new_start_line, 0,
                    msg='Start line out of bounds for ' + where)
                self.assertLessEqual(
                    new_end_line, slice1.size,
                    msg='End line out of bounds for ' + where)
                slice2 = slice1[new_start_line:new_end_line + 1]

                # Check results
                self.assertTrue(np.all(ref == slice2),
                                msg='Incorrect slice for ' + where)

    def test_save_gac(self):
        """Test selection of user defined scanlines

        Scanline Nr:  1  2  3  4  5  6  7  8  9 10 11 12 13 14
        Missing:      X           X                 X        X
        Inv. lat/lon:    X  X  X                 X     X  X

        => All missing lines: [1, 2, 3, 4, 5, 10, 11, 12, 13, 14
        """
        import pygac.gac_io

        # Define test data
        fv = pygac.gac_io.MISSING_DATA_LATLON
        along = 10
        across = 3
        miss_lines = np.array([1, 5, 11, 14])  # never recorded
        qualflags = np.array(
            [[2, 3, 4, 6, 7, 8, 9, 10, 12, 13]]).transpose()  # scanline number
        lats = np.arange(along*across).reshape((along, across))
        lats[:2+1, :] = fv  # Make first 3 lines invalid
        lats[7:, :] = fv  # Make last 3 lines invalid
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
            lats_ref = lats[new_start_line:new_end_line + 1, :]
            qualflags_ref = qualflags[new_start_line:new_end_line+1]
            xutcs_ref = xutcs[new_start_line:new_end_line + 1]
            midnight_scanline_ref = midnight_scanlines_ref[itest]

            with mock.patch('pygac.gac_io.avhrrGAC_io') as io_mock:
                pygac.gac_io.save_gac(
                    satellite_name='dummy',
                    xutcs=xutcs,
                    lats=lats,
                    lons=lons,
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
                    mask=np.zeros(lats.shape, dtype='bool'),
                    qual_flags=qualflags,
                    start_line=start_line,
                    end_line=end_line,
                    tsmcorr=False,
                    gac_file='dummy',
                    midnight_scanline=midnight_scanline,
                    miss_lines=miss_lines,
                    switch=np.zeros(lats.shape, dtype='bool')
                )
                expected_args = [
                    mock.ANY,
                    CalledWithAdapter(xutcs_ref),
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    CalledWithAdapter(lats_ref*1000.0),
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
                    CalledWithAdapter(qualflags_ref),
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    mock.ANY,
                    midnight_scanline_ref,
                    CalledWithAdapter(all_miss_lines_ref)
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
