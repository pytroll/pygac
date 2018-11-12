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
from pygac.gac_io import update_start_end_line


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
        import pygac.gac_io

        # Define test data
        fv = pygac.gac_io.MISSING_DATA_LATLON
        along = 10
        across = 3
        lats = np.arange(along*across).reshape((along, across))
        lats[:2+1, :] = fv  # Make first 3 lines invalid
        lats[7:, :] = fv  # Make last 3 lines invalid
        lons = lats.copy()
        xutcs = np.arange(along).astype('datetime64[ms]')
        qualflags = np.array(
            [[2, 3, 4, 6, 7, 8, 9, 10, 12, 13]]).transpose()  # scanline number
        miss_lines = np.array([1, 5, 11, 14])
        midnight_scanline = 5
        dummydata = np.ones(lats.shape)

        for start_line in (0, 1, 4):
            for end_line in (9, 8, 5):
                where = 'with start/end line {0}/{1}'.format(start_line,
                                                             end_line)

                # Define reference data: Ignore lines with invalid lat/lon info
                # (here 0, 1, 2, 7, 8, 9)
                new_start_line = max(3, start_line)
                new_end_line = min(6, end_line)
                lats_ref = lats[new_start_line:new_end_line + 1, :]
                qualflags_ref = qualflags[new_start_line:new_end_line+1]
                qualflags_pre = qualflags
                xutcs_ref = xutcs[new_start_line:new_end_line + 1]

                # Midnight scanline has to be updated to new slice
                midnight_scanline_ref = midnight_scanline - new_start_line

                # Define Tester
                def test(satellite_name, xutcs, startdate, enddate, starttime,
                         endtime, lats, lons, ref1, ref2, ref3, bt3, bt4, bt5,
                         sun_zen, sat_zen, sun_azi, sat_azi, rel_azi,
                         qual_flags, start_line, end_line,
                         total_number_of_scan_lines, last_scan_line_number,
                         corr, gac_file, midnight_scanline, miss_lines):
                    """Check whether scanlines with invalid lat/lon info have
                     been excluded correctly"""

                    # Rescale lats
                    lats[np.where(lats != fv)] /= 1000

                    # Test whether the arrays have been sliced correctly
                    self.assertTrue(
                        np.all(lats == lats_ref),
                        msg='Latitude has not been sliced correctly ' + where
                    )
                    self.assertTrue(
                        np.all(qual_flags == qualflags_ref),
                        msg='Qualflags have not been sliced correctly ' + where
                    )
                    self.assertTrue(
                        np.all(xutcs == xutcs_ref),
                        msg='UTC times have not been sliced correctly ' + where
                    )

                    # If missing lines are tracked correctly, the union of
                    # pre-slicing scanline numbers and post-slicing missing
                    # lines should include all scanlines
                    self.assertEqual(
                        set(qualflags_pre[:, 0]).union(set(miss_lines)),
                        set(range(1, 14+1)),
                        msg='Missing scanlines have not been tracked '
                            'correctly ' + where
                    )

                    # Check whether midnigt scanline has been updated correctly
                    self.assertEqual(
                        midnight_scanline, midnight_scanline_ref,
                        msg='Midnight scanline has not been updated '
                            'correctly ' + where
                    )

                # Patch gac_io.avhrrGAC_io using the above tester
                pygac.gac_io.avhrrGAC_io = test

                # Call save_gac which then calls the tester instead of
                # gac_io.avhrrGAC_io
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

        # Delete the module as we don't want the patched methods to affect
        # other tests
        del pygac.gac_io


def suite():
    """The suite for test_io"""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestIO))
    return mysuite


if __name__ == '__main__':
    unittest.main()
