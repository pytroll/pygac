#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013, 2014 Martin Raspaud

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
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

"""Test the GAC KLM reader."""

import datetime as dt
from unittest import mock

import numpy as np
import numpy.testing

from pygac.gac_klm import GACKLMReader
from pygac.klm_reader import header
from pygac.lac_klm import LACKLMReader, scanline
from pygac.tests.utils import CalledWithArray


class TestKLM:
    """Test the klm reader."""

    def setup(self):
        """Set up the tests."""
        self.reader = GACKLMReader()

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
        assert time == dt.datetime(2019, 5, 3, 0, 2, 3, 456000)

    def test_get_times(self):
        """Test readout of scanline timestamps."""
        self.reader.scans = {'scan_line_year': 1,
                             'scan_line_day_of_year': 2,
                             'scan_line_utc_time_of_day': 3}
        assert self.reader._get_times() == (1, 2, 3)

    def test_get_ch3_switch(self):
        """Test channel 3 identification."""
        self.reader.scans = {
            'scan_line_bit_field': np.array([1, 2, 3, 4, 5, 6])}
        switch_exp = np.array([1, 2, 3, 0, 1, 2])
        numpy.testing.assert_array_equal(
            self.reader.get_ch3_switch(), switch_exp)

    def test_postproc(self):
        """Test KLM specific postprocessing."""
        self.reader.scans = {
            'scan_line_bit_field': np.array([0, 1, 2])}
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

    def test_quality_indicators(self):
        """Test the quality indicator unpacking."""
        reader = self.reader
        QFlag = reader.QFlag
        quality_indicators = np.array([
            0,  # nothing flagged
            -1,  # everything flagged
            QFlag.CALIBRATION | QFlag.NO_EARTH_LOCATION,
            QFlag.TIME_ERROR | QFlag.DATA_GAP,
            QFlag.FATAL_FLAG
        ], dtype=np.uint32)
        reader.scans = {self.reader._quality_indicators_key: quality_indicators}
        # test mask, i.e. QFlag.FATAL_FLAG | QFlag.CALIBRATION | QFlag.NO_EARTH_LOCATION
        expected_mask = np.array([False, True, True, False, True], dtype=bool)
        numpy.testing.assert_array_equal(reader.mask, expected_mask)
        # test individual flags
        assert reader._get_corrupt_mask(flags=QFlag.FATAL_FLAG).any()
        # count the occurence (everything flagged and last entrance => 2)
        assert reader._get_corrupt_mask(flags=QFlag.FATAL_FLAG).sum() == 2


class TestGACKLM:
    """Tests for gac klm."""

    def setup(self):
        """Set up the tests."""
        self.reader = GACKLMReader()

    @mock.patch('pygac.klm_reader.get_tsm_idx')
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


class TestLACKLM:
    """Tests for lac klm."""

    def setup(self):
        """Set up the tests."""
        self.reader = LACKLMReader()
        self.reader.scans = np.ones(100, dtype=scanline)
        self.reader.head = np.ones(1, dtype=header)[0]
        self.reader.spacecraft_id = 12
        self.reader.head["noaa_spacecraft_identification_code"] = self.reader.spacecraft_id
        self.reader.spacecraft_name = "metopa"
        self.reader.scans["scan_line_number"] = np.arange(100)
        # PRT
        self.reader.scans["telemetry"]["PRT"] = 400
        self.reader.scans["telemetry"]["PRT"][0::5, :] = 0

    def test_get_ch3_switch(self):
        """Test channel 3 identification."""
        self.reader.scans = {
            'scan_line_bit_field': np.array([1, 2, 3, 4, 5, 6])}
        switch_exp = np.array([1, 2, 3, 0, 1, 2])
        numpy.testing.assert_array_equal(
            self.reader.get_ch3_switch(), switch_exp)

    def test_calibrate_channels(self):
        """Test channel calibration."""
        # ICT
        self.reader.scans["back_scan"] = 400
        self.reader.scans["back_scan"][0::5, :] = 0
        # Space
        self.reader.scans["space_data"] = 400
        self.reader.scans["space_data"][0::5, :] = 0

        assert np.any(np.isfinite(self.reader.get_calibrated_channels()))

    def test_calibrate_inactive_3b(self):
        """Test calibration of an inactive 3b."""
        channels = self.reader.get_calibrated_channels()
        assert np.all(np.isnan(channels[:, :, 3]))
