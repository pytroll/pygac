#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013, 2014 Martin Raspaud

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>

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

"""Test the GAC KLM reader.
"""

import datetime as dt
import numpy as np
import numpy.testing
import unittest
import subprocess as sp
try:
    from unittest import mock
except ImportError:
    import mock

from pygac import gac_klm
from pygac.tests.utils import CalledWithArray


# test_file = "test_data/NSS.GHRR.NL.D01183.S1313.E1458.B0400708.GC"
# ref_result = "test_data/ECC_avhrrGAC_noaa16_20010702Z131338_20010702Z145809.h5"
# my_result = "ECC_avhrrGAC_noaa16_20010702Z131338_20010702Z145809.h5"

# ref_sunsat = "test_data/ECC_sunsatGAC_noaa16_20010702Z131338_20010702Z145809.h5.orig"
# my_sunsat = "ECC_sunsatGAC_noaa16_20010702Z131338_20010702Z145809.h5"

test_file = "test_data/NSS.GHRR.NL.D02187.S1904.E2058.B0921517.GC"
ref_result = "test_data/ECC_GAC_sunsatangles_noaa16_99999_20020706T1904020Z_20020706T2058095Z.h5"
my_result = "/tmp/ECC_GAC_sunsatangles_noaa16_99999_20020706T1904020Z_20020706T2058095Z.h5"

ref_sunsat = "test_data/ECC_GAC_sunsatangles_noaa16_99999_20020706T1904020Z_20020706T2058095Z.h5"
my_sunsat = "/tmp/ECC_GAC_sunsatangles_noaa16_99999_20020706T1904020Z_20020706T2058095Z.h5"


class TestKLM(unittest.TestCase):
    def setUp(self):
        self.reader = gac_klm.GACKLMReader()

    def test_global(self):
        gac_klm.main(test_file, 0, 0)

        child = sp.Popen(["h5diff", ref_sunsat, my_sunsat],
                         stdout=sp.PIPE)
        streamdata = child.communicate()[0]
        retc = child.returncode
        self.assertTrue(retc == 0, msg=streamdata)

        child = sp.Popen(["h5diff", ref_result, my_result],
                         stdout=sp.PIPE)
        streamdata = child.communicate()[0]
        retc = child.returncode
        self.assertTrue(retc == 0, msg=streamdata)

    def test_get_lonlat(self):
        """Test readout of lon/lat coordinates."""
        earth_loc = 1e4 * np.array([[1, 2, 3, 4],
                                    [5, 6, 7, 8]])
        self.reader.scans = {'earth_location': earth_loc}

        lons_exp = np.array([[2, 4],
                             [6, 8]])
        lats_exp = np.array([[1, 3],
                             [5, 7]])

        lons, lats = self.reader._get_lonlat()
        numpy.testing.assert_array_equal(lons, lons_exp)
        numpy.testing.assert_array_equal(lats, lats_exp)

    def test_get_header_timestamp(self):
        """Test readout of header timestamp."""
        self.reader.head = {
            'start_of_data_set_year': np.array([2019]),
            'start_of_data_set_day_of_year': np.array([123]),
            'start_of_data_set_utc_time_of_day': np.array([123456])
        }
        time = self.reader.get_header_timestamp()
        self.assertEqual(time, dt.datetime(2019, 5, 3, 0, 2, 3, 456000))

    def test_get_times(self):
        """Test readout of scanline timestamps."""
        self.reader.scans = {'scan_line_year': 1,
                             'scan_line_day_of_year': 2,
                             'scan_line_utc_time_of_day': 3}
        self.assertTupleEqual(self.reader._get_times(), (1, 2, 3))

    def test_get_ch3_switch(self):
        """Test channel 3 identification."""
        self.reader.scans = {
            'scan_line_bit_field': np.array([1, 2, 3, 4, 5, 6])}
        switch_exp = np.array([1, 2, 3, 0, 1, 2])
        numpy.testing.assert_array_equal(
            self.reader.get_ch3_switch(), switch_exp)

    @mock.patch('pygac.gac_klm.GACKLMReader.get_ch3_switch')
    def test_postproc(self, get_ch3_switch):
        """Test KLM specific postprocessing."""
        get_ch3_switch.return_value = np.array([0, 1, 2])
        channels = np.array([[[1., 2., 3., 4.],
                              [1., 2., 3., 4.]],
                             [[1., 2., 3., 4.],
                              [1., 2., 3., 4.]],
                             [[1., 2., 3, 4.],
                              [1., 2., 3, 4.]]])  # (lines, pixels, channels)

        masked_exp = np.array([[[1., 2., np.nan, 4.],
                                [1., 2., np.nan, 4.]],
                               [[1., 2., 3., np.nan],
                                [1., 2., 3., np.nan]],
                               [[1., 2., np.nan, np.nan],
                                [1., 2., np.nan, np.nan]]])
        self.reader.postproc(channels)  # masks in-place
        numpy.testing.assert_array_equal(channels, masked_exp)

    @mock.patch('pygac.gac_klm.get_tsm_idx')
    def test_get_tsm_pixels(self, get_tsm_idx):
        """Test channel set used for TSM correction."""
        ones = np.ones((409, 100))
        zeros = np.zeros(ones.shape)
        ch1 = 1*ones
        ch2 = 2*ones
        ch4 = 4*ones
        ch5 = 5*ones
        channels = np.dstack((ch1, ch2, zeros, zeros, ch4, ch5))
        self.reader.get_tsm_pixels(channels)
        get_tsm_idx.assert_called_with(CalledWithArray(ch1),
                                       CalledWithArray(ch2),
                                       CalledWithArray(ch4),
                                       CalledWithArray(ch5))


if __name__ == '__main__':
    unittest.main()
