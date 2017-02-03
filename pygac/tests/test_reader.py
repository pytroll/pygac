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
import numpy as np
from pygac.gac_reader import GACReader


class GACReaderDummy(GACReader):
    def _get_times(self):
        pass

    def get_header_timestamp(self):
        pass

    def read(self, filename):
        pass

    @property
    def tsm_affected_intervals(self):
        pass


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

    def test_midnight_scanline(self):
        """Test midnight scanline computation"""
        reader = GACReaderDummy()

        # Define test cases...
        # ... midnight scanline exists
        utcs1 = (24 * 3600 * 1000 + np.array([-3, -2, -1, 0, 1, 2, 3])).astype(
            'datetime64[ms]')
        scanline1 = 2

        # ... midnight scanline does not exist
        utcs2 = np.array([1, 2, 3]).astype('datetime64[ms]')
        scanline2 = None

        for utcs, scanline in zip((utcs1, utcs2), (scanline1, scanline2)):
            reader.utcs = utcs
            reader.times = utcs.astype(datetime.datetime)
            self.assertEqual(reader.get_midnight_scanline(), scanline,
                             msg='Incorrect midnight scanline')


if __name__ == '__main__':
    unittest.main()
