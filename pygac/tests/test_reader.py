#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author(s):

#   Stephan Finkensieper <stephan.finkensieper@dwd.de>
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
"""Test the readers."""

import datetime
import logging
import os
import sys
import unittest
from unittest import mock

import numpy as np
import numpy.testing
import pytest
import xarray as xr

from pygac.gac_pod import scanline
from pygac.gac_pod import scanline as gacpod_scanline
from pygac.gac_reader import GACReader, ReaderError
from pygac.lac_pod import LACPODReader
from pygac.lac_pod import scanline as lacpod_scanline
from pygac.lac_reader import LACReader
from pygac.pod_reader import POD_QualityIndicator, header3
from pygac.pod_reader import tbm_header as tbm_header_dtype
from pygac.reader import (
    MAX_SCAN_ANGLES,
    NAVIGATION_METADATA_SCHEMA_VERSION,
    NoTLEData,
    frame_buffer_lag_for,
    has_a_trusted_clock,
    max_scan_angle_for,
    seconds_from_scanlines,
    should_fit_the_clock,
    should_fit_the_pitch,
    yaw_steers,
)


class FakePath(os.PathLike):
    """Fake path class."""

    def __init__(self, path):
        """Initialize the path."""
        self.path = str(path)

    def __fspath__(self):
        """Return the path."""
        return self.path


class FakeGACReader_withtimes(GACReader):
    """Fake GAC reader class."""

    QFlag = POD_QualityIndicator
    _quality_indicators_key = "quality_indicators"
    tsm_affected_intervals = {None: []}
    along_track = 3
    across_track = 409

    def __init__(self, *args, **kwargs):
        """Initialize the fake reader."""
        super().__init__(*args, **kwargs)
        self.scan_width = self.across_track
        scans = np.zeros(self.along_track, dtype=scanline)
        scans["scan_line_number"] = np.arange(self.along_track)
        scans["sensor_data"] = 128
        self.scans = scans
        self.head = {"foo": "bar"}
        self.head = np.rec.fromrecords([("bar", ),],names="foo")
        self.spacecraft_name = "noaa14"
        self.spacecrafts_orbital = {3: "noaa14"}
        self.spacecraft_id = 3

    def _get_times_from_file(self):
        # 2000-11-17T03:18:44.806
        year = np.full(self.along_track, 2000, dtype=int)
        jday = np.full(self.along_track, 322, dtype=int)
        msec = 1000 * 11924.806 * np.arange(1, self.along_track+1, dtype=int)
        return year, jday, msec

    def get_header_timestamp(self):
        """Get the header timestamp."""
        return datetime.datetime(2000, 11, 17, 3, 18, 44, 803)

    def get_telemetry(self):
        """Get the telemetry.

        Contract, matching PODReader/KLMReader.get_telemetry:
            mean_prt_counts   (scan_line_index,)
            full_space_counts (scan_line_index, 10, n_channels)
            full_ict_counts   (scan_line_index, 10, 3)
        """
        prt = 51 * np.ones(self.along_track)  # prt threshold is 50
        prt[::5] = 0
        full_space = 101 * np.ones((self.along_track, 10, 5))  # space threshold is 100
        full_ict = 101 * np.ones((self.along_track, 10, 3))  # ict threshold is 100
        return prt, full_space, full_ict

    def get_vis_telemetry(self):
        """Get the telemetry."""
        space = 40 * np.ones((self.along_track, 3))  # space threshold is 100
        total_space = 40 * np.ones((self.along_track, 10, 3))  # space threshold is 100
        return space, total_space

    def _adjust_clock_drift(self):
        pass

    def _get_lonlat_from_file(self):
        return np.zeros((self.along_track, 51)), np.zeros((self.along_track, 51))

    @staticmethod
    def _get_ir_channels_to_calibrate():
        return [3, 4, 5]

    def postproc(self, channels):
        """Postprocess the data."""
        pass

    def read(self, filename, fileobj=None):
        """Read the data."""
        pass

    @classmethod
    def read_header(cls, filename, fileobj=None):
        """Read the header."""
        pass

    def get_tsm_pixels(self, channels):
        """Get the tsm pixels."""
        pass

    #def get_tle_lines(self):
    #    """Get tle lines."""
    #    pass

class FakeGACReader(GACReader):
    """Fake GAC reader class."""

    QFlag = POD_QualityIndicator
    _quality_indicators_key = "quality_indicators"
    tsm_affected_intervals = {None: []}
    along_track = 3
    across_track = 409

    def __init__(self, *args, **kwargs):
        """Initialize the fake reader."""
        super().__init__(*args, **kwargs)
        self.scan_width = self.across_track
        scans = np.zeros(self.along_track, dtype=scanline)
        scans["scan_line_number"] = np.arange(self.along_track)
        scans["sensor_data"] = 128
        self.scans = scans
        self.head = {"foo": "bar"}
        self.head = np.rec.fromrecords([("bar", ),],names="foo")
        self.spacecraft_name = "noaa6"

    def _get_times_from_file(self):
        year = np.full(self.along_track, 1970, dtype=int)
        jday = np.full(self.along_track, 1, dtype=int)
        msec = 1000 * np.arange(1, self.along_track+1, dtype=int)
        return year, jday, msec

    def get_header_timestamp(self):
        """Get the header timestamp."""
        return datetime.datetime(1970, 1, 1)

    def get_telemetry(self):
        """Get the telemetry.

        Contract, matching PODReader/KLMReader.get_telemetry:
            mean_prt_counts   (scan_line_index,)
            full_space_counts (scan_line_index, 10, n_channels)
            full_ict_counts   (scan_line_index, 10, 3)
        """
        prt = 51 * np.ones(self.along_track)  # prt threshold is 50
        prt[::5] = 0
        full_space = 101 * np.ones((self.along_track, 10, 5))  # space threshold is 100
        full_ict = 101 * np.ones((self.along_track, 10, 3))  # ict threshold is 100
        return prt, full_space, full_ict

    def get_vis_telemetry(self):
        """Get the telemetry."""
        space = 40 * np.ones((self.along_track, 3))  # space threshold is 100
        total_space = 40 * np.ones((self.along_track, 10, 3))  # space threshold is 100
        return space, total_space

    def _adjust_clock_drift(self):
        pass

    def _get_lonlat_from_file(self):
        return np.zeros((self.along_track, 51)), np.zeros((self.along_track, 51))

    @staticmethod
    def _get_ir_channels_to_calibrate():
        return [3, 4, 5]

    def postproc(self, channels):
        """Postprocess the data."""
        pass

    def read(self, filename, fileobj=None):
        """Read the data."""
        pass

    @classmethod
    def read_header(cls, filename, fileobj=None):
        """Read the header."""
        pass

    def get_tsm_pixels(self, channels):
        """Get the tsm pixels."""
        pass

    def get_tle_lines(self):
        """Get tle lines."""
        pass

class FakeGACReaderWithWrongPRTs(FakeGACReader_withtimes):

    along_track = 10

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.scans["scan_line_number"] = np.array([1, 3, 4, 5, 6, 7, 8, 9, 10, 11])

@pytest.fixture
def pod_file_with_tbm_header(tmp_path):
    """Create a pod file with a tbm header, and header, and some scanlines."""
    number_of_scans = 3

    tbm_header = np.zeros(1, dtype=tbm_header_dtype)
    tbm_header["data_set_name"] = "BRN.HRPT.NJ.D00322.S0334.E0319.B3031919.BL\x80\x80".encode("cp500")
    tbm_header["select_flag"] = b"S"
    tbm_header["beginning_latitude"] = b"+77"
    tbm_header["ending_latitude"] = b"+22"
    tbm_header["beginning_longitude"] = b"-004"
    tbm_header["ending_longitude"] = b"+032"
    tbm_header["start_hour"] = b"AL"
    tbm_header["start_minute"] = b"L "
    tbm_header["number_of_minutes"] = b"ALL"
    tbm_header["appended_data_flag"] = b"Y"
    tbm_header["channel_select_flag"][0, :5] = b"\x01"
    tbm_header["sensor_data_word_size"] = b"10"

    header = np.zeros(1, dtype=header3)
    header["noaa_spacecraft_identification_code"] = 3
    header["data_type_code"] = 48
    header["start_time"] = [51522, 181, 62790]
    header["number_of_scans"] = number_of_scans
    header["end_time"] = [51522, 195, 42286]
    header["processing_block_id"] = b"3031919"
    header["ramp_auto_calibration"] = 0
    header["number_of_data_gaps"] = 0
    header["dacs_quality"] = [21, 7, 0, 0, 0, 0]
    header["calibration_parameter_id"] = 12336
    header["dacs_status"] = 24
    header["nadir_earth_location_tolerance"] = 20
    header["start_of_data_set_year"] = 2000
    # EBCDIC, aka cp500 encoding
    header["data_set_name"] = (b"\xc2\xd9\xd5K\xc8\xd9\xd7\xe3K\xd5\xd1K\xc4\xf0\xf0\xf3\xf2\xf2K\xe2\xf0\xf3\xf1\xf9K"
                               b"\xc5\xf0\xf3\xf3\xf4K\xc2\xf3\xf0\xf3\xf1\xf9\xf1\xf9K\xc2\xd3@@")
    header["year_of_epoch"] = 2000
    header["julian_day_of_epoch"] = 322
    header["millisecond_utc_epoch_time_of_day"] = 11864806
    header["semi_major_axis"] = 7226239
    header["eccentricity"] = 100496
    header["inclination"] = 9915900
    header["argument_of_perigee"] = 2511798
    header["right_ascension"] = 30366227
    header["mean_anomaly"] = 7341437
    header["x_component_of_position_vector"] = 40028760
    header["y_component_of_position_vector"] = -60151283
    header["z_component_of_position_vector"] = 0,
    header["x_dot_component_of_position_vector"] = -990271
    header["y_dot_component_of_position_vector"] = -646066
    header["z_dot_component_of_position_vector"] = 7337861
    header["yaw_fixed_error_correction"] = 0
    header["roll_fixed_error_correction"] = 0
    header["pitch_fixed_error_correction"] = 0

    scanlines = np.zeros(number_of_scans, dtype=lacpod_scanline)
    scanlines["scan_line_number"] = np.arange(number_of_scans) + 1
    scanlines["time_code"][:, :2] = [51522, 181]
    scanlines["time_code"][:, 2] = np.arange(62790, 62790 + 166 * number_of_scans, 166)
    scanlines["quality_indicators"] = 1073741824
    scanlines["calibration_coefficients"] = [149722896, -23983736, 173651248, -27815402, -1803673, 6985664, -176321840,
                                             666680320, -196992576, 751880640]
    scanlines["number_of_meaningful_zenith_angles_and_earth_location_appended"] = 51

    one_line_angles = np.array([-24, -26, -28, -30, -32, -33, -34, -35, -37, -37, -38, -39, -40, -41, -41, -42,
                                -43, -43, -44, -45, -45, -46, -46, -47, -48, -48, -49, -49, -50, -50, -51, -52,
                                -52, -53, -53, -54, -55, -55, -56, -57, -58, -59, -60, -61, -62, -63, -64, -66,
                                -68, -70, -73])
    scanlines["solar_zenith_angles"] = [one_line_angles, one_line_angles + 1, one_line_angles + 2]

    one_line_lats = np.array([9929, 9966, 9983, 9988, 9984, 9975, 9963, 9948, 9931, 9913, 9895, 9875,
                              9856, 9836, 9816, 9796, 9775, 9755, 9734, 9713, 9692, 9671, 9650, 9628,
                              9606, 9583, 9560, 9537, 9513, 9488, 9463, 9437, 9409, 9381, 9351, 9320,
                              9288, 9253, 9216, 9177, 9135, 9090, 9040, 8986, 8926, 8859, 8784, 8697,
                              8596, 8476, 8327])

    scanlines["earth_location"]["lats"] = [one_line_lats - 1, one_line_lats, one_line_lats + 1]

    one_line_lons = np.array([-36, 747, 1415, 1991, 2492, 2931, 3320, 3667, 3979, 4262, 4519, 4755, 4972,
                              5174, 5362, 5538, 5703, 5860, 6009, 6150, 6286, 6416, 6542, 6663, 6781,
                              6896, 7009, 7119, 7227, 7334, 7440, 7545, 7650, 7755, 7860, 7966, 8073,
                              8181, 8291, 8403, 8518, 8637, 8759, 8886, 9019, 9159, 9307, 9466, 9637,
                              9824, 10033])

    scanlines["earth_location"]["lons"] = [one_line_lons, one_line_lons + 1, one_line_lons + 2]

    scanlines["telemetry"] = 2047
    scanlines["sensor_data"] = 99
    scanlines["add_on_zenith"] = 0
    scanlines["clock_drift_delta"] = 0
    scanlines["sensor_data"] = 2047

    pod_filename = tmp_path / "image.l1b"
    offset = 14800
    fill = np.zeros(offset - header3.itemsize, dtype=np.uint8)

    with open(pod_filename, "wb") as fd_:
        fd_.write(tbm_header.tobytes())
        fd_.write(header.tobytes())
        fd_.write(fill.tobytes())
        fd_.write(scanlines.tobytes())
    return pod_filename


@pytest.fixture
def pod_file_with_tbm_header_gac(tmp_path):
    """Create a pod file (GAC) with a tbm header, and header, and some scanlines."""
    number_of_scans = 3

    tbm_header = np.zeros(1, dtype=tbm_header_dtype)
    tbm_header["data_set_name"] = "BRN.GHRR.NJ.D00322.S0334.E0319.B3031919.BL\x80\x80".encode("cp500")
    tbm_header["select_flag"] = b"S"
    tbm_header["beginning_latitude"] = b"+77"
    tbm_header["ending_latitude"] = b"+22"
    tbm_header["beginning_longitude"] = b"-004"
    tbm_header["ending_longitude"] = b"+032"
    tbm_header["start_hour"] = b"AL"
    tbm_header["start_minute"] = b"L "
    tbm_header["number_of_minutes"] = b"ALL"
    tbm_header["appended_data_flag"] = b"Y"
    tbm_header["channel_select_flag"][0, :5] = b"\x01"
    tbm_header["sensor_data_word_size"] = b"10"

    header = np.zeros(1, dtype=header3)
    header["noaa_spacecraft_identification_code"] = 3
    header["data_type_code"] = 32
    header["start_time"] = [51522, 181, 62790]
    header["number_of_scans"] = number_of_scans
    header["end_time"] = [51522, 195, 42286]
    header["processing_block_id"] = b"3031919"
    header["ramp_auto_calibration"] = 0
    header["number_of_data_gaps"] = 0
    header["dacs_quality"] = [21, 7, 0, 0, 0, 0]
    header["calibration_parameter_id"] = 12336
    header["dacs_status"] = 24
    header["nadir_earth_location_tolerance"] = 20
    header["start_of_data_set_year"] = 2000
    # EBCDIC, aka cp500 encoding
    header["data_set_name"] = (b"\xc2\xd9\xd5K\xc8\xd9\xd7\xe3K\xd5\xd1K\xc4\xf0\xf0\xf3\xf2\xf2K\xe2\xf0\xf3\xf1\xf9K"
                               b"\xc5\xf0\xf3\xf3\xf4K\xc2\xf3\xf0\xf3\xf1\xf9\xf1\xf9K\xc2\xd3@@")
    header["year_of_epoch"] = 2000
    header["julian_day_of_epoch"] = 322
    header["millisecond_utc_epoch_time_of_day"] = 11864806
    header["semi_major_axis"] = 7226239
    header["eccentricity"] = 100496
    header["inclination"] = 9915900
    header["argument_of_perigee"] = 2511798
    header["right_ascension"] = 30366227
    header["mean_anomaly"] = 7341437
    header["x_component_of_position_vector"] = 40028760
    header["y_component_of_position_vector"] = -60151283
    header["z_component_of_position_vector"] = 0,
    header["x_dot_component_of_position_vector"] = -990271
    header["y_dot_component_of_position_vector"] = -646066
    header["z_dot_component_of_position_vector"] = 7337861
    header["yaw_fixed_error_correction"] = 0
    header["roll_fixed_error_correction"] = 0
    header["pitch_fixed_error_correction"] = 0

    scanlines = np.zeros(number_of_scans, dtype=gacpod_scanline)
    scanlines["scan_line_number"] = np.arange(number_of_scans) + 1
    scanlines["time_code"][:, :2] = [51522, 181]
    scanlines["time_code"][:, 2] = np.arange(62790, 62790 + 166 * number_of_scans, 166)
    scanlines["quality_indicators"] = 1073741824
    scanlines["calibration_coefficients"] = [149722896, -23983736, 173651248, -27815402, -1803673, 6985664, -176321840,
                                             666680320, -196992576, 751880640]
    scanlines["number_of_meaningful_zenith_angles_and_earth_location_appended"] = 51

    one_line_angles = np.array([-24, -26, -28, -30, -32, -33, -34, -35, -37, -37, -38, -39, -40, -41, -41, -42,
                                -43, -43, -44, -45, -45, -46, -46, -47, -48, -48, -49, -49, -50, -50, -51, -52,
                                -52, -53, -53, -54, -55, -55, -56, -57, -58, -59, -60, -61, -62, -63, -64, -66,
                                -68, -70, -73])
    scanlines["solar_zenith_angles"] = [one_line_angles, one_line_angles + 1, one_line_angles + 2]

    one_line_lats = np.array([9929, 9966, 9983, 9988, 9984, 9975, 9963, 9948, 9931, 9913, 9895, 9875,
                              9856, 9836, 9816, 9796, 9775, 9755, 9734, 9713, 9692, 9671, 9650, 9628,
                              9606, 9583, 9560, 9537, 9513, 9488, 9463, 9437, 9409, 9381, 9351, 9320,
                              9288, 9253, 9216, 9177, 9135, 9090, 9040, 8986, 8926, 8859, 8784, 8697,
                              8596, 8476, 8327])

    scanlines["earth_location"]["lats"] = [one_line_lats - 1, one_line_lats, one_line_lats + 1]

    one_line_lons = np.array([-36, 747, 1415, 1991, 2492, 2931, 3320, 3667, 3979, 4262, 4519, 4755, 4972,
                              5174, 5362, 5538, 5703, 5860, 6009, 6150, 6286, 6416, 6542, 6663, 6781,
                              6896, 7009, 7119, 7227, 7334, 7440, 7545, 7650, 7755, 7860, 7966, 8073,
                              8181, 8291, 8403, 8518, 8637, 8759, 8886, 9019, 9159, 9307, 9466, 9637,
                              9824, 10033])

    scanlines["earth_location"]["lons"] = [one_line_lons, one_line_lons + 1, one_line_lons + 2]

    scanlines["telemetry"] = 2047
    scanlines["sensor_data"] = 99
    scanlines["add_on_zenith"] = 0
    scanlines["clock_drift_delta"] = 0
    scanlines["sensor_data"] = 2047

    pod_filename = tmp_path / "image.l1b"
    offset = 14800
    fill = np.zeros(offset - header3.itemsize, dtype=np.uint8)

    with open(pod_filename, "wb") as fd_:
        fd_.write(tbm_header.tobytes())
        fd_.write(header.tobytes())
        fd_.write(fill.tobytes())
        fd_.write(scanlines.tobytes())
    return pod_filename


@pytest.fixture()
def pod_tle(tmp_path):
    lines = ("1 23455U 94089A   00322.04713399  .00000318  00000-0  19705-3 0  5298\n"
             "2 23455  99.1591 303.5706 0010037  25.7760 334.3905 14.12496755303183\n"
             "1 23455U 94089A   00322.96799836  .00000229  00000-0  14918-3 0  5303\n"
             "2 23455  99.1590 304.5117 0009979  23.1101 337.0518 14.12496633303313\n")
    tle_filename = tmp_path / "noaa14.tle"
    with tle_filename.open("w") as fd:
        fd.write(lines)
    return tle_filename


def test_get_calibrated_channels(pod_file_with_tbm_header_gac,pod_tle):
    """Test getting calibrated channels."""
    reader = FakeGACReader_withtimes(pod_file_with_tbm_header_gac,
                                     tle_dir=pod_tle.parent,
                                     tle_name=pod_tle.name)
    reader.interpolate_coords = True
    res = reader.get_calibrated_channels()

    # Old test results
    #np.testing.assert_allclose(res[:, 1, 0], 8.84714652)
    #np.testing.assert_allclose(res[:, 2, 1], 10.23511303)
    np.testing.assert_allclose(res[:, 1, 0], 11.24028967)
    np.testing.assert_allclose(res[:, 2, 1], 13.98060798)
    assert reader.meta_data["calib_coeffs_version"] == "PATMOS-x, v2023"


def test_get_calibrated_channels_with_wrong_prts(pod_file_with_tbm_header_gac,
                                                 pod_tle):
    """Test getting calibrated channels."""
    reader = FakeGACReaderWithWrongPRTs(pod_file_with_tbm_header,
                                        tle_dir=pod_tle.parent,
                                        tle_name=pod_tle.name)
    reader.interpolate_coords = True
    res = reader.get_calibrated_channels()

    # Old test results
    #np.testing.assert_allclose(res[:, 1, 0], 8.84714652)
    #np.testing.assert_allclose(res[:, 2, 1], 10.23511303)
    np.testing.assert_allclose(res[:, 1, 0], 11.24028967)
    np.testing.assert_allclose(res[:, 2, 1], 13.98060798)
    assert reader.meta_data["calib_coeffs_version"] == "PATMOS-x, v2023"


def preset_times(reader, utcs):
    """Install fixed scanline timestamps on *reader*.

    Since pygac commit 0d73911 ``Reader.get_times`` is ``@cache``-decorated and no
    longer short-circuits when ``_times_as_np_datetime64`` is already set, so it
    unconditionally re-reads and re-corrects the timestamps from the file. A test
    that wants specific timestamps must therefore also stop it doing that.
    """
    reader._times_as_np_datetime64 = utcs
    reader.get_times = lambda: reader._times_as_np_datetime64
    return utcs


class TestGacReader(unittest.TestCase):
    """Test the common GAC Reader."""

    longMessage = True

    @mock.patch.multiple("pygac.gac_reader.GACReader",
                         __abstractmethods__=set())
    def setUp(self, *mocks):
        """Set up the tests."""
        self.reader = GACReader()

    def test_filename(self):
        """Test the setter of the filename property."""
        # test path with .gz extension
        filename = "NSS.GHRR.TN.D80001.S0332.E0526.B0627173.WI"
        filepath = "/path/to/" + filename + ".gz"
        self.reader.filename = filepath
        self.assertEqual(self.reader.filename, filename)
        self.reader.filename = None
        self.assertIsNone(self.reader.filename)
        self.reader.filename = FakePath(filepath)
        self.assertEqual(self.reader.filename, filename)

    @unittest.skipIf(sys.version_info.major < 3, "Skipped in python2!")
    def test__read_scanlines(self):
        """Test the scanline extraction."""
        self.reader.scanline_type = np.dtype([
            ("a", "S2"), ("b", "<i4")])
        # request more scan lines than available
        with self.assertWarnsRegex(RuntimeWarning,
                                   "Unexpected number of scanlines!"):
            buffer = (b"a\x00\x01\x00\x00\x00"
                      b"b\x00\x02\x00\x00\x00"
                      b"c\x00\x03\x00\x00\x00")
            count = 4
            self.reader._read_scanlines(buffer, count)
        # check the output
        first_line = self.reader.scans[0]
        self.assertEqual(first_line["a"], b"a")
        self.assertEqual(first_line["b"], 1)

    def test__validate_header(self):
        """Test the header validation."""
        # wrong name pattern
        with self.assertRaisesRegex(ReaderError,
                                    "Data set name .* does not match!"):
            head = {"data_set_name": b"abc.txt"}
            self.reader._validate_header(head)
        # Unicode errors are now caught with the same exception.
        with self.assertRaisesRegex(ReaderError,
                                    "Data set name .* does not match!"):
            head = {"data_set_name": b"\xea\xf8"}
            self.reader._validate_header(head)

    def test__correct_data_set_name_ebcdic_encoded_header_invalid_path(self):
        """Test the data_set_name correction in file header."""
        inv_filename = "InvalidFileName"
        inv_filepath = "path/to/" + inv_filename

        expected_data_set_name = "NSS.GHRR.TN.D80001.S0332.E0526.B0627173.WI"
        val_head = {"data_set_name": "NSS.GHRR.TN.D80001.S0332.E0526.B0627173.WI".encode("cp500")}
        head = self.reader._correct_data_set_name(val_head.copy(), inv_filepath)
        assert head["data_set_name"] == expected_data_set_name.encode()

    def test__correct_data_set_name_valid_header_and_file(self):
        """Test the data_set_name correction in file header."""
        val_filename = "NSS.GHRR.TN.D80001.S0332.E0526.B0627173.WI"
        val_filepath = "path/to/" + val_filename
        val_head = {"data_set_name": b"NSS.GHRR.TN.D80001.S0332.E0526.B0627173.WI"}
        # Note: always pass a copy to _correct_data_set_name, because
        #       the input header is modified in place.
        # enter a valid data_set_name and filepath
        head = self.reader._correct_data_set_name(val_head.copy(), val_filepath)
        assert head["data_set_name"] == val_filename.encode()

    def test__correct_data_set_name_invalid_header_and_valid_file(self):
        """Test the data_set_name correction in file header."""
        val_filename = "NSS.GHRR.TN.D80001.S0332.E0526.B0627173.WI"
        val_filepath = "path/to/" + val_filename
        inv_head = {"data_set_name": b"InvalidDataSetName"}

        # enter an invalid data_set_name, but valid filepath
        head = self.reader._correct_data_set_name(inv_head.copy(), val_filepath)
        assert head["data_set_name"] == val_filename.encode()

    def test__correct_data_set_name_invalid_header_and_file(self):
        """Test the data_set_name correction in file header."""
        inv_filename = "InvalidFileName"
        inv_filepath = "path/to/" + inv_filename
        inv_head = {"data_set_name": b"InvalidDataSetName"}
        with self.assertRaisesRegex(ReaderError, "Cannot determine data_set_name!"):
            _ = self.reader._correct_data_set_name(inv_head.copy(), inv_filepath)

    def test__correct_data_set_name_valid_header_invalid_file(self):
        """Test the data_set_name correction in file header."""
        val_head = {"data_set_name": b"NSS.GHRR.TN.D80001.S0332.E0526.B0627173.WI"}
        inv_filename = "InvalidFileName"
        inv_filepath = "path/to/" + inv_filename

        # enter a valid data_set_name, and an invalid filepath
        # should be fine, because the data_set_name is the pefered source
        head = self.reader._correct_data_set_name(val_head.copy(), inv_filepath)
        assert head["data_set_name"] == val_head["data_set_name"]

    def test__correct_data_set_name_valid_header_pathlib_file(self):
        """Test the data_set_name correction in file header."""
        val_filename = "NSS.GHRR.TN.D80001.S0332.E0526.B0627173.WI"
        val_filepath = "path/to/" + val_filename
        val_head = {"data_set_name": b"NSS.GHRR.TN.D80001.S0332.E0526.B0627173.WI"}

        fs_filepath = FakePath(val_filepath)
        head = self.reader._correct_data_set_name(val_head.copy(), fs_filepath)
        self.assertEqual(head["data_set_name"], val_filename.encode())

    @mock.patch("pygac.reader.Reader.get_calibrated_channels")
    def test__get_calibrated_channels_uniform_shape(self, get_channels):
        """Test the uniform shape as required by gac_io.save_gac."""
        # check if it raises the assertion error
        channels = np.arange(2*2*5, dtype=float).reshape((2, 2, 5))
        get_channels.return_value = channels
        with self.assertRaises(AssertionError):
            self.reader._get_calibrated_channels_uniform_shape()

    def test_to_datetime64(self):
        """Test conversion from (year, jday, msec) to datetime64."""
        t0 = GACReader.to_datetime64(year=np.array(1970), jday=np.array(1),
                                     msec=np.array(0))
        self.assertEqual(t0.astype("i8"), 0,
                         msg="Conversion (year, jday, msec) to datetime64 "
                             "is not correct")

    def test_to_datetime(self):
        """Test conversion from datetime64 to datetime."""
        dt = datetime.datetime(2019, 10, 1, 12, 15, 18, 12345)
        dt64 = np.datetime64(dt)
        self.assertEqual(self.reader.to_datetime(dt64), dt)

    def test_lineno2msec(self):
        """Test scanline timestamp estimation."""
        self.assertEqual(self.reader.lineno2msec(12345), 6172000)

    @mock.patch("pygac.reader.Reader.QFlag", POD_QualityIndicator)
    @mock.patch("pygac.reader.Reader.update_meta_data")
    @mock.patch("pygac.gac_reader.GACReader._get_lonlat_from_file")
    @mock.patch("pygac.gac_reader.GACReader._get_corrupt_mask")
    @mock.patch("pygac.gac_reader.GACReader._adjust_clock_drift")
    @mock.patch("pygac.reader.Reader.lonlat_interpolator")
    def test_get_lonlat(self, lonlat_interpolator, adjust_clockdrift,
                        get_corrupt_mask, get_lonlat,
                        update_meta_data):
        """Test common lon/lat computation."""
        lon_i = np.array([np.nan, 1, 2, 3, -180.1, 180.1])
        lat_i = np.array([1, 2, 3, np.nan, -90.1, 90.1])
        get_lonlat.return_value = lon_i, lat_i
        lonlat_interpolator.return_value = lon_i, lat_i
        get_corrupt_mask.return_value = np.array(
            [0, 0, 1, 0, 0, 0], dtype=bool)

        lons_exp = np.array([np.nan, 1, np.nan, 3., np.nan, np.nan])
        lats_exp = np.array([1, 2, np.nan, np.nan, np.nan, np.nan])

        # Default
        lons, lats = self.reader.get_lonlat()
        get_lonlat.assert_called()
        update_meta_data.assert_called()
        adjust_clockdrift.assert_called()
        numpy.testing.assert_array_equal(lons, lons_exp)
        numpy.testing.assert_array_equal(lats, lats_exp)

        # Interpolation disabled
        lonlat_interpolator.reset_mock()
        adjust_clockdrift.reset_mock()
        self.reader.clock_drift_correction_applied = False
        self.reader.interpolate_coords = False
        self.reader.adjust_clock_drift = True
        self.reader.lons = self.reader.lats = None
        self.reader.get_lonlat()
        numpy.testing.assert_array_equal(lons, lons_exp)
        numpy.testing.assert_array_equal(lats, lats_exp)
        lonlat_interpolator.assert_not_called()
        adjust_clockdrift.assert_called()

        # Clock drift adjustment disabled
        lonlat_interpolator.reset_mock()
        adjust_clockdrift.reset_mock()
        self.reader.interpolate_coords = True
        self.reader.adjust_clock_drift = False
        self.reader.lons = self.reader.lats = None
        self.reader.get_lonlat()
        numpy.testing.assert_array_equal(lons, lons_exp)
        numpy.testing.assert_array_equal(lats, lats_exp)
        lonlat_interpolator.assert_called()
        adjust_clockdrift.assert_not_called()

        # Test caching
        methods = [get_lonlat, lonlat_interpolator,
                   adjust_clockdrift, get_corrupt_mask]
        for method in methods:
            method.reset_mock()
        self.reader.get_lonlat()
        for method in methods:
            method.asser_not_called()

    @mock.patch("pygac.reader.Reader.QFlag", POD_QualityIndicator)
    @mock.patch("pygac.reader.Reader.update_meta_data")
    @mock.patch("pygac.gac_reader.GACReader._get_corrupt_mask")
    @mock.patch("pygac.gac_reader.GACReader._adjust_clock_drift")
    @mock.patch("pygac.gac_reader.GACReader._get_lonlat_from_file")
    @mock.patch("pygac.reader.Reader.lonlat_interpolator")
    def test_interpolate(self, lonlat_interpolator, _get_lonlat, _adjust_clock_drift,
                         _get_corrupt_mask, update_meta_data):
        """Test interpolate method in get_lonlat."""
        self.lons = None
        self.lats = None
        rng = np.random.default_rng()
        lr_lons = 90 * rng.random((17, 51))
        lr_lats = 90 * rng.random((17, 51))
        _get_lonlat.return_value = lr_lons, lr_lats
        self.interpolate_coords = True
        lonlat_interpolator.reset_mock()
        lonlat_interpolator.return_value = (90 * rng.random((17, 409)),
                                          90 * rng.random((17, 409)))
        lons, lats = self.reader.get_lonlat()
        self.assertEqual(lons.shape[1], 409)
        lonlat_interpolator.assert_called_once_with(lr_lons, lr_lats)

    @mock.patch("pygac.gac_reader.GACReader._get_corrupt_mask")
    def test_get_corrupt_mask(self, get_corrupt_mask):
        """Test common computation of corrupt scanline mask."""
        get_corrupt_mask.return_value = [1, 2, 3]
        self.assertEqual(self.reader.mask, [1, 2, 3])

        # Test caching
        get_corrupt_mask.reset_mock()
        self.reader.mask
        get_corrupt_mask.assert_not_called()

    def test_midnight_scanline(self):
        """Test midnight scanline computation."""
        # Define test cases...
        # ... midnight scanline exists
        utcs1 = np.array([-3, -2, -1, 0, 1, 2, 3]).astype("datetime64[ms]")
        scanline1 = 2

        # ... midnight scanline does not exist
        utcs2 = np.array([1, 2, 3]).astype("datetime64[ms]")
        scanline2 = None

        for utcs, scan_line in zip((utcs1, utcs2), (scanline1, scanline2)):
            preset_times(self.reader, utcs)
            self.assertEqual(self.reader.get_midnight_scanline(), scan_line,
                             msg="Incorrect midnight scanline")

    def test_miss_lines(self):
        """Test detection of missing scanlines."""
        lines = [2, 4, 5, 6, 10, 11, 12]
        miss_lines_ref = [1, 3, 7, 8, 9]
        self.reader.scans = np.zeros(
            len(lines), dtype=[("scan_line_number", "i2")])
        self.reader.scans["scan_line_number"] = lines
        miss_lines = self.reader.get_miss_lines()
        self.assertTrue((miss_lines == miss_lines_ref).all(),
                        msg="Missing scanlines not detected correctly")
        self.assertEqual(miss_lines.dtype, int)

    def test_tle2datetime64(self, *mocks):
        """Test conversion from TLE timestamps to datetime64."""
        dates = np.array([70365.1234, 18001.25])
        dates64_exp = [np.datetime64(datetime.datetime(1970, 12, 31, 2, 57, 41, 760000), "[ms]"),
                       np.datetime64(datetime.datetime(2018, 1, 1, 6, 0), "[ms]")]
        dates64 = GACReader.tle2datetime64(dates)
        self.assertTrue(np.all(dates64 == dates64_exp))

    @mock.patch("pygac.gac_reader.GACReader.get_times")
    @mock.patch("pygac.gac_reader.GACReader.get_tle_file")
    @mock.patch("pygac.gac_reader.GACReader.read_tle_file")
    def test_get_tle_lines(self, read_tle_file, *mocks):
        """Test identification of closest TLE lines."""
        tle_data = ["1 38771U 12049A   18363.63219793 -.00000013  00000-0  14176-4 0  9991\r\n",
                    "2 38771  98.7297  60.1350 0002062  95.9284  25.0713 14.21477560325906\r\n",
                    "1 38771U 12049A   18364.62426010 -.00000015  00000-0  13136-4 0  9990\r\n",  # 2018-12-30 14:58
                    "2 38771  98.7295  61.1159 0002062  94.5796  60.2561 14.21477636326047\r\n",
                    "1 38771U 12049A   18364.94649306 -.00000018  00000-0  12040-4 0  9996\r\n",  # 2018-12-30 22:42
                    "2 38771  98.7295  61.4345 0002060  94.1226 268.7521 14.21477633326092\r\n",
                    "1 38771U 12049A   18365.81382142 -.00000015  00000-0  13273-4 0  9991\r\n",
                    "2 38771  98.7294  62.2921 0002057  92.7653  26.0030 14.21477711326215\r\n"]

        expected = {
            datetime.datetime(2018, 12, 20, 12, 0): None,
            datetime.datetime(2018, 12, 28, 12, 0): 0,
            datetime.datetime(2018, 12, 30, 16, 0): 2,
            datetime.datetime(2018, 12, 30, 20, 0): 4,
            datetime.datetime(2019, 1, 1, 12, 0): 6,
            datetime.datetime(2019, 1, 8, 12, 0): None
        }

        read_tle_file.return_value = tle_data
        for time, tle_idx in expected.items():
            self.reader._times_as_np_datetime64 = np.array([time], dtype="datetime64[ms]")
            self.reader.tle_lines = None
            if tle_idx is None:
                self.assertRaises(NoTLEData, self.reader.get_tle_lines)
            else:
                tle1, tle2 = self.reader.get_tle_lines()
                self.assertEqual(tle1, tle_data[tle_idx])
                self.assertEqual(tle2, tle_data[tle_idx + 1])

    @mock.patch("pygac.reader.Reader.get_tle_lines")
    def test_get_sat_angles_without_tle_at_nadir(self, get_tle_lines):
        """Test that the get satellite angles without tle at nadir."""
        get_tle_lines.side_effect = NoTLEData("No TLE data available")
        rng = np.random.RandomState(125)
        self.reader.lons = rng.rand(100, 409) * 90
        self.reader.lats = rng.rand(100, 409) * 90
        self.reader._times_as_np_datetime64 = np.array(
            [numpy.datetime64(datetime.datetime(1980, 1, 3, 11, 47, 15, 469000)) for date in range(100)])
        sat_azi, sat_elev = self.reader.get_sat_angles()
        self.assertEqual(np.sum(np.isnan(sat_elev)), 0)
        np.testing.assert_allclose(sat_elev[:, 204], 90., atol=0.01)

    @mock.patch("pygac.reader.Reader.get_tle_lines")
    def test_get_sat_angles_without_tle(self, get_tle_lines):
        """Test the get satellite angles without tle."""
        get_tle_lines.side_effect = NoTLEData("No TLE data available")

        # Test data correspond to columns 0:2, 201:208 and 407:409. Extracted like this:
        # self.lons[0:5, [0, 1, 201, 202, 203, 204, 205, 206, 207, -2, -1]]
        self.reader.lons = np.array([[69.41555135, 68.76815744, 28.04133742, 27.94671757, 27.85220562,
                                      27.7578125, 27.66354783, 27.5694182, 27.47542957, 2.66416611,
                                      2.39739436],
                                     [69.41409536, 68.76979616, 28.00228658, 27.9076628, 27.8131467,
                                      27.71875, 27.62448295, 27.53035209, 27.43636312, 2.61727408,
                                      2.35049275],
                                     [69.42987929, 68.78543423, 27.96407251, 27.86936406, 27.77457923,
                                      27.6796875, 27.58467527, 27.48959853, 27.39453053, 2.5704025,
                                      2.30362323],
                                     [69.44430772, 68.80064104, 27.91910034, 27.82340242, 27.72794715,
                                      27.6328125, 27.53805662, 27.44366144, 27.34959008, 2.53088093,
                                      2.26729486],
                                     [69.47408815, 68.8259859, 27.87666513, 27.78224611, 27.68795045,
                                      27.59375, 27.49962326, 27.40557682, 27.31162435, 2.48359319,
                                      2.21976689]])
        self.reader.lats = np.array([[71.62830288, 71.67081539, 69.90976034, 69.89297223, 69.87616536,
                                      69.859375, 69.84262997, 69.82593315, 69.80928089, 61.61334632,
                                      61.466097],
                                     [71.65222644, 71.69455292, 69.9331981, 69.91640991, 69.89960295,
                                      69.8828125, 69.86606741, 69.84937054, 69.83271828, 61.62903602,
                                      61.48182119],
                                     [71.68324681, 71.7256586, 69.96444365, 69.94765637, 69.93085095,
                                      69.9140625, 69.89731955, 69.8806245, 69.86397337, 61.65247459,
                                      61.50526043],
                                     [71.71489701, 71.75715746, 69.98789144, 69.97110185, 69.9542928,
                                      69.9375, 69.92075278, 69.90405431, 69.88740114, 61.66489162,
                                      61.51597801],
                                     [71.7390608, 71.78108857, 70.01126173, 69.99449392, 69.97770882,
                                      69.9609375, 69.9442054, 69.92751555, 69.91086544, 61.67769951,
                                      61.52802445]])
        self.reader._times_as_np_datetime64 = np.array(
            [numpy.datetime64(datetime.datetime(1980, 1, 3, 11, 47, 15, 469000)),
             numpy.datetime64(datetime.datetime(1980, 1, 3, 11, 47, 15, 969000)),
             numpy.datetime64(datetime.datetime(1980, 1, 3, 11, 47, 16, 469000)),
             numpy.datetime64(datetime.datetime(1980, 1, 3, 11, 47, 16, 969000)),
             numpy.datetime64(datetime.datetime(1980, 1, 3, 11, 47, 17, 469000))])
        expected_sat_azi_0 = np.array([283.09872924, 283.12775589, 283.13951497, 283.14786413, 283.19638805])
        expected_sat_azi_201 = np.array([272.85051989, 273.79847634, 272.04794616, 273.01363377, 274.00055397])
        expected_sat_azi_408 = np.array([39.77021472, 39.71516966, 39.68104134, 39.60503726, 39.5403431])
        expected_sat_elev_0 = np.array([20.94889204, 20.96041284, 20.96368521, 20.96826309, 20.95849212])
        expected_sat_elev_204 = np.array([89.24533884, 89.22663677, 89.25079817, 89.24938043, 89.23004118])
        sat_azi, sat_elev = self.reader.get_sat_angles()
        np.testing.assert_allclose(sat_azi[:, 0], expected_sat_azi_0, atol=1.0)
        np.testing.assert_allclose(sat_azi[:, 2], expected_sat_azi_201, atol=35.0)  # Azi bad close to center!
        np.testing.assert_allclose(sat_azi[:, -1], expected_sat_azi_408, atol=1.0)
        np.testing.assert_allclose(sat_elev[:, 0], expected_sat_elev_0, atol=1.0)
        np.testing.assert_allclose(sat_elev[:, 5], expected_sat_elev_204, atol=1.0)

    @mock.patch("pygac.gac_reader.GACReader._get_corrupt_mask")
    def test_get_angles(self, get_corrupt_mask):
        """Test get_angles function of the reader."""
        # Line: 1, 649, 6198 and 12658 from Tiros-N file (1980-01-03 11:47)
        lon_i = np.array(
            [69.41555, 152.10587, 164.3131, 67.23855, np.nan])[:, np.newaxis]
        lat_i = np.array(
            [71.6283, 85.24265, -62.076958, 82.72296, np.nan])[:, np.newaxis]
        get_corrupt_mask.return_value = np.isnan(lon_i)
        self.reader.lons = lon_i
        self.reader.lats = lat_i
        self.reader.tle_lines = [
            '1 11060U 78096A   80003.54792075  .00000937  00000-0  52481-3 0  2588\r\n',  # noqa
            '2 11060  98.9783 332.1605 0012789  88.8047 271.4583 14.11682873 63073\r\n']  # noqa
        preset_times(self.reader, np.array(
            [315748035469, 315748359969,
             315751135469, 315754371969,
             315754371969]).astype("datetime64[ms]"))
        self.reader.spacecrafts_orbital = {25: "tiros n"}
        self.reader.spacecraft_id = 25
        expected_sat_azi = np.array(
            [-76.90, 11.08, 145.33, -50.01, np.nan])[:, np.newaxis]
        expected_sun_azi = np.array(
            [-120.36, -31.94, -173.51, -93.67, np.nan])[:, np.newaxis]
        expected_sat_zenith = np.array(
            [69.05, 69.04, 69.55, 69.07, np.nan])[:, np.newaxis]
        expected_sun_zenith = np.array(
            [104.30, 116.94, 94.86, 112.60, np.nan])[:, np.newaxis]
        expected_rel_azi = np.array(
            [43.45, 43.01, 41.16, 43.65, np.nan])[:, np.newaxis]

        retv = self.reader.get_angles()
        (sat_azi, sat_zenith, sun_azi, sun_zenith, rel_azi) = retv
        np.testing.assert_allclose(sat_azi, expected_sat_azi, atol=0.01)
        np.testing.assert_allclose(sun_azi, expected_sun_azi, atol=0.01)
        np.testing.assert_allclose(sat_zenith, expected_sat_zenith, atol=0.01)
        np.testing.assert_allclose(sun_zenith, expected_sun_zenith, atol=0.01)
        np.testing.assert_allclose(rel_azi, expected_rel_azi, atol=0.01)

    def test_get_tle_file(self):
        """Test get_tle_file."""
        self.reader.tle_dir = "/tle/dir"
        self.reader.tle_name = "tle_%(satname)s.txt"
        self.reader.spacecraft_name = "ISS"
        tle_file = self.reader.get_tle_file()
        self.assertEqual(tle_file, "/tle/dir/tle_ISS.txt")

    @mock.patch("pygac.gac_reader.GACReader.get_tsm_pixels")
    def test_mask_tsm_pixels(self, get_tsm_pixels):
        """Test masking of pixels affected by the scan motor issue."""
        get_tsm_pixels.return_value = ([0, 1], [0, 1])
        channels = np.array([[[1., 2., 3.],
                              [4., 5., 6.]],
                             [[7., 8., 9.],
                              [10., 11., 12.]]])  # (lines, pixels, channels)
        masked_exp = np.array([[[np.nan, np.nan, np.nan],
                                [4., 5., 6.]],
                               [[7., 8., 9.],
                                [np.nan, np.nan, np.nan]]])

        channels = xr.DataArray(channels, dims=["scan_line_index", "columns", "channel_name"],
                                coords=dict(channel_name=["1", "2", "3a"] ))
        ds = xr.Dataset(dict(channels=channels))

        self.reader.mask_tsm_pixels(ds)  # masks in-place
        numpy.testing.assert_array_equal(ds["channels"].values, masked_exp)

    def test_correct_scan_line_numbers(self):
        """Test scanline number correction."""
        scans, expected = _get_scanline_numbers(14000)
        self.reader.scans = scans
        self.reader.correct_scan_line_numbers()
        numpy.testing.assert_array_equal(self.reader.scans["scan_line_number"],
                                         expected)

    @mock.patch("pygac.gac_reader.GACReader.get_header_timestamp")
    def test_correct_times_thresh(self, get_header_timestamp):
        """Test correction of scanline timestamps."""
        header_time = datetime.datetime(2016, 8, 16, 16, 7, 36)

        # Create artificial timestamps
        _, scan_line_numbers = _get_scanline_numbers(14000)
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
        self.reader._times_as_np_datetime64 = msecs.astype(">M8[ms]")
        self.reader.scans = np.array(scan_line_numbers,
                                     dtype=[("scan_line_number", ">u2")])

        # Test correction
        self.reader.correct_times_thresh()
        numpy.testing.assert_array_equal(self.reader._times_as_np_datetime64, utcs_expected)

    def test_calculate_sun_earth_distance_correction(self):
        """Test the calculate sun earth distance correction method."""
        preset_times(self.reader, np.array([315748035469, 315748359969,
                                     315751135469, 315754371969,
                                     315754371969]).astype("datetime64[ms]"))
        corr = self.reader.get_sun_earth_distance_correction()
        numpy.testing.assert_almost_equal(corr, 0.96660494, decimal=7)

    @mock.patch("pygac.reader.Reader.get_sun_earth_distance_correction")
    @mock.patch("pygac.reader.Reader.get_midnight_scanline")
    @mock.patch("pygac.reader.Reader.get_miss_lines")
    def test_update_metadata(self,
                             get_miss_lines,
                             get_midnight_scanline,
                             get_sun_earth_distance_correction):
        """Test updating the metadata."""
        get_miss_lines.return_value = "miss_lines"
        get_midnight_scanline.return_value = "midn_line"
        get_sun_earth_distance_correction.return_value = "factor"
        self.reader.head = {"foo": "bar"}

        self.reader.update_meta_data()

        mda_exp = {"midnight_scanline": "midn_line",
                   "missing_scanlines": "miss_lines",
                   "sun_earth_distance_correction_factor": "factor",
                   "gac_header": {"foo": "bar"}}
        self.assertDictEqual(self.reader.meta_data, mda_exp)


def _get_scanline_numbers(scanlines_along_track):
    """Create artificial scanline numbers with some corruptions.

    Returns:
        Corrupted and corrected scanline numbers.

    """
    scans = np.zeros(scanlines_along_track, dtype=[("scan_line_number", ">u2")])
    scans["scan_line_number"] = np.arange(1, scanlines_along_track + 1)

    # ... with 500 missing scanlines at scanline 8000
    scans["scan_line_number"][8000:] += 500
    corrected = scans["scan_line_number"].copy()

    # ... and some spikes here and there
    scans["scan_line_number"][3000] += 1E4
    scans["scan_line_number"][9000] -= 1E4
    corrected = np.delete(corrected, [3000, 9000])

    return scans, corrected


class TestLacReader(unittest.TestCase):
    """Test the common LAC Reader."""

    longMessage = True

    @mock.patch.multiple("pygac.lac_reader.LACReader",
                         __abstractmethods__=set())
    def setUp(self, *mocks):
        """Set up the tests."""
        self.reader = LACReader()

    def test_lac_reader_accepts_FRAC(self):
        """Test the header validation."""
        head = {"data_set_name": b"NSS.FRAC.M1.D19115.S2352.E0050.B3425758.SV"}
        self.reader._validate_header(head)

    def test_correct_scan_line_numbers(self):
        """Test scanline number correction."""
        scans, expected = _get_scanline_numbers(22000)
        self.reader.scans = scans
        self.reader.correct_scan_line_numbers()
        numpy.testing.assert_array_equal(self.reader.scans["scan_line_number"],
                                         expected)


def test_podlac_eosip(pod_file_with_tbm_header):
    """Test reading a real podlac file."""
    reader = LACPODReader(interpolate_coords=False)
    reader.read(pod_file_with_tbm_header)
    assert reader.head.itemsize == header3.itemsize
    # this is broken in eosip pod data, tbm data set name has start and end times reversed.
    # assert reader.head["data_set_name"] == reader.tbm_head["data_set_name"]
    # todo: test that duplicate lines are removed and recomputed


def test_read_to_dataset_is_a_dataset_including_channels_and_telemetry(pod_file_with_tbm_header, pod_tle):
    """Test creating an xr.Dataset from a gac file."""
    import xarray as xr
    reader = LACPODReader(interpolate_coords=True, tle_dir=pod_tle.parent, tle_name=pod_tle.name)
    dataset = reader.read_as_dataset(pod_file_with_tbm_header)
    assert isinstance(dataset, xr.Dataset)
    assert dataset["channels"].shape == (3, 2048, 5)
    assert dataset.attrs["processing_block_id"] == b"3031919"
    assert "times" in dataset.coords
    assert "scan_line_index" in dataset.coords
    assert "channel_name" in dataset.coords
    assert dataset["mean_prt_counts"].shape == (3,)
    assert dataset["full_ict_counts"].shape == (3, 10, 3)
    assert dataset["full_space_counts"].shape == (3, 10, 5)
    assert dataset["quality_flags"].shape == (3, )
    assert "flag_meanings" in dataset["quality_flags"].attrs
    assert "flag_masks" in dataset["quality_flags"].attrs
#
# Code seems to require interpolation to with the new addition
# of get_angles in getting the counts, so this test doesn't work
# Comment by J.Mittaz, UoR
#
#def test_read_to_dataset_without_interpolation(pod_file_with_tbm_header, pod_tle):
#    """Test creating an xr.Dataset from a gac file."""
#    reader = LACPODReader(interpolate_coords=False, tle_dir=pod_tle.parent, tle_name=pod_tle.name)
#    dataset = reader.read_as_dataset(pod_file_with_tbm_header)
#
#    assert dataset["longitude"].shape == (3, 51)
#    assert dataset["latitude"].shape == (3, 51)

def test_read_to_dataset_with_interpolation(pod_file_with_tbm_header,
                                            pod_tle):
    """Test creating an xr.Dataset from a gac file."""
    reader = LACPODReader(interpolate_coords=True, tle_dir=pod_tle.parent, tle_name=pod_tle.name)
    dataset = reader.read_as_dataset(pod_file_with_tbm_header)

    assert dataset.coords["longitude"].shape == (3, 2048)
    assert dataset.coords["latitude"].shape == (3, 2048)


def test_passing_calibration_coeffs_to_reader_init_is_deprecated():
    """Test that passing calibration-specific info to the reader's init is not allowed anymore."""
    with pytest.warns(PendingDeprecationWarning):
        FakeGACReader(calibration_method="noaa", custom_calibration=dict(a="1"))
    with pytest.warns(PendingDeprecationWarning):
        FakeGACReader(calibration_method="noaa", calibration_file="somefile")


def test_passing_calibration_to_reader(pod_file_with_tbm_header,pod_tle):
    """Test passing calibration info to `get_calibrated_channels`."""
    method = "InvalidMethod"
    with pytest.raises(ValueError, match=method):
        reader = FakeGACReader_withtimes(
                               calibration_method=method,
                               tle_dir=pod_tle.parent,
                               tle_name=pod_tle.name,
                               interpolate_coords=True)

    reader = FakeGACReader_withtimes(calibration_method="noaa",
                           tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                           interpolate_coords=True)

    res = reader.get_calibrated_channels()

    # np.testing.assert_allclose(res[:, 1, 0], 8.84714652)
    # np.testing.assert_allclose(res[:, 2, 1], 10.23511303)
    np.testing.assert_allclose(res[:, 1, 0], 11.24028967)
    np.testing.assert_allclose(res[:, 2, 1], 13.98060798)
    assert reader.meta_data["calib_coeffs_version"] == "PATMOS-x, v2023"


def test_computing_lonlats(pod_file_with_tbm_header, pod_tle):
    """Test computing lons and lats from TLE data."""
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True)
    dataset = reader.read_as_dataset(pod_file_with_tbm_header)
    lons = dataset["longitude"].values
    assert lons[0, 0] == pytest.approx(-4.5185248)
    assert lons[0, -1] == pytest.approx(79.5159661)


def test_recomputing_lonlats_with_time_offset(pod_file_with_tbm_header, pod_tle):
    """Test computing lons and lats from TLE data with time offset."""
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name, compute_lonlats_from_tles=True)
    reader.read(pod_file_with_tbm_header)
    lons, lats = reader._compute_lonlats(time_offset=np.timedelta64(500, "ms"))
    assert lons[0, 0] == pytest.approx(-4.4772295)
    assert lons[0, -1] == pytest.approx(79.4561871)


def test_georeferencing_with_first_guess(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """Test getting a first guess for georeferencing from the file lon/lats."""
    expected_time_offset = 2  # seconds
    def skip_thermal(channels, *args, **kwargs):
        # This fixture has only three scan lines and no usable PRT readings,
        # so the thermal calibration cannot run on it; these tests are about
        # georeferencing.
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    def mock_disp(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        return 0, (0, 0, 0), ([10000] * 60, [1000] * 60)
    use_a_stub_georeferencer(monkeypatch, mock_disp)
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name, compute_lonlats_from_tles=True,
                          reference_image="some_world_image.tif", adjust_clock_drift=False)
    def file_lonlats(self, *args, **kwargs):
        lons, lats = self._compute_lonlats(time_offset=np.timedelta64(int(expected_time_offset * 1e9), "ns"))
        return lons[:, self.lonlat_sample_points], lats[:, self.lonlat_sample_points]
    reader._get_lonlat_from_file = file_lonlats.__get__(reader)
    reader.read(pod_file_with_tbm_header)
    start_time = reader.get_times()[0]
    _ = reader.get_calibrated_dataset()
    new_start_time = reader.get_times()[0]
    tdiff = start_time - new_start_time + np.timedelta64(expected_time_offset, "s")
    # The pre-alignment fits one time offset to four corners, so how exactly it
    # recovers a planted one follows the geolocation underneath it. This fixture
    # carries three scan lines, where that fit is at its least constrained.
    assert tdiff < np.timedelta64(50, "ms")


def test_georeferencing_fails(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """Test georeferencing."""

    def skip_thermal(channels, *args, **kwargs):
        # This fixture has only three scan lines and no usable PRT readings,
        # so the thermal calibration cannot run on it; these tests are about
        # georeferencing.
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    def mock_disp(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        return 0, (0, 0, 0), ([10000], [10000])
    use_a_stub_georeferencer(monkeypatch, mock_disp)
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name, compute_lonlats_from_tles=True,
                          reference_image="some_world_image.tif")
    reader.read(pod_file_with_tbm_header)
    with pytest.warns(RuntimeWarning):
        _ = reader.get_calibrated_dataset()


def use_a_stub_georeferencer(monkeypatch, stub):
    """Install a stub for the two calls pygac makes, given one written as a single answer.

    pygac measures the displacement, decides from it whether the platform's clock is worth
    fitting, and only then fits. A test that only cares what comes back from the fit says so
    once, here, rather than spelling both halves out.
    """
    from georeferencer import georeferencer

    def measures(calibrated_ds, *args, **rest):
        measures.answer = stub(calibrated_ds, *args, **rest)
        # as many control points as the stub recorded, as the real georeferencer returns
        count = calibrated_ds["gcp_x"].size
        return np.column_stack([np.linspace(0.0, 2.0, count), np.linspace(60.0, 1990.0, count)]), \
            np.zeros((count, 2)), 0.0

    def fits(gcps, ref_lons, ref_lats, start_time, tle, max_scan_angle, **rest):
        return measures.answer

    from pyorbital import geoloc_avhrr
    monkeypatch.setattr(georeferencer, "measure_swath_displacement", measures)
    monkeypatch.setattr(geoloc_avhrr, "fit_navigation", answering_as_a_navigation_fit(fits))


def answering_as_a_navigation_fit(fits):
    """Wrap a stub written as the old fit's answer so it answers as pyorbital's fit does.

    The stubs say what came back -- time, attitude, and the distances before and after --
    which is all most of these tests care about. The fit now also hands back what it left
    of each point's miss and how far it can be trusted; a stub that says nothing of those
    leaves the misses at zero and the uncertainty unknown.
    """
    from pyorbital.geoloc_avhrr import NavigationFit

    def fit_navigation(*args, **kwargs):
        time_offset, attitude, (_, after) = fits(*args, **kwargs)
        after = np.asarray(after, dtype=float)
        # one miss per control point, as pyorbital gives them; the count the stub reports
        # is carried by the distances, which is what the reader counts
        return NavigationFit(time_offset, tuple(attitude), np.zeros((len(args[0]), 2)), after,
                             np.full((4, 4), np.nan), {})

    return fit_navigation

def navigate_with(monkeypatch, count=60, time_offset=0.0, attitude=(0.0, 0.0, 0.0), distance=1000.0,
                  misses=None, covariance=None, settings=None, along_track_lines=0.0, gcp_lonlats=None,
                  fits=None, field=None):
    """Stand in for the matching and for pyorbital's fit, saying everything the fit returns.

    The matching leaves a coherent field behind and reports *along_track_lines* of coarse
    shift; the fit leaves every point *distance* metres off, with *misses* in swath pixels
    and the parameters' *covariance* in seconds and radians.
    """
    from georeferencer import georeferencer
    from pyorbital import geoloc_avhrr
    from pyorbital.geoloc_avhrr import NavigationFit
    fits_asked = fits

    def measures(calibrated_ds, *args, **rest):
        (field or record_a_coherent_field)(calibrated_ds, count)
        spread_over_the_swath = np.column_stack([np.linspace(0.0, 2.0, count), np.linspace(60.0, 1990.0, count)])
        lonlats = np.zeros((count, 2)) if gcp_lonlats is None else gcp_lonlats
        return spread_over_the_swath, lonlats, along_track_lines

    def fits(gcps, *args, **rest):
        if fits_asked is not None:
            fits_asked.append(rest)
        return NavigationFit(time_offset, tuple(attitude),
                             np.zeros((count, 2)) if misses is None else misses,
                             np.full(count, float(distance)),
                             np.zeros((4, 4)) if covariance is None else covariance,
                             {"loss": "soft_l1"} if settings is None else settings)

    monkeypatch.setattr(georeferencer, "measure_swath_displacement", measures)
    monkeypatch.setattr(geoloc_avhrr, "fit_navigation", fits)

def record_a_coherent_field(calibrated_ds, count=60):
    """Populate the control point arrays the way the real georeferencer does.

    Its displacements are what the acceptance checks read, so a double standing in
    for it has to leave them behind too.
    """
    spread = np.linspace(0.0, 2000.0, count)
    for name, values in (("gcp_y", spread), ("gcp_x", spread),
                         ("gcp_y_displacement", np.full(count, 3.0)),
                         ("gcp_x_displacement", np.full(count, -2.0))):
        calibrated_ds[name] = xr.DataArray(values, dims=["points"])


def test_georeferencing(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """Test georeferencing."""

    def skip_thermal(channels, *args, **kwargs):
        # This fixture has only three scan lines and no usable PRT readings,
        # so the thermal calibration cannot run on it; these tests are about
        # georeferencing.
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    def mock_disp(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        return 0.5, (0, 0, 0), ([10000] * 60, [1000] * 60)
    use_a_stub_georeferencer(monkeypatch, mock_disp)
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name, compute_lonlats_from_tles=True,
                          reference_image="some_world_image.tif")
    reader.read(pod_file_with_tbm_header)
    dataset = reader.get_calibrated_dataset()
    assert dataset.attrs["max_scan_angle"] == 55.37
    lons = dataset["longitude"].values
    assert lons[0, 0] == pytest.approx(-4.4772295)
    assert lons[0, -1] == pytest.approx(79.4561871)

def test_orthocorrection(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """Test computing lons and lats from TLE data."""

    def skip_thermal(channels, *args, **kwargs):
        # This fixture has only three scan lines and no usable PRT readings,
        # so the thermal calibration cannot run on it; these tests are about
        # georeferencing.
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    def mock_disp(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        return 0.5, (0, 0, 0), ([10000] * 60, [1000] * 60)
    use_a_stub_georeferencer(monkeypatch, mock_disp)
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name, compute_lonlats_from_tles=True,
                          reference_image="some_world_image.tif", dem="dem_file.tif")
    reader.read(pod_file_with_tbm_header)


    def mock_orthocorrection(calibrated_ds, *args):
        calibrated_ds["tc_lons"] = calibrated_ds["longitude"] + 0.0001
        calibrated_ds["tc_lats"] = calibrated_ds["latitude"] + 0.0001
        return calibrated_ds

    from georeferencer import georeferencer
    monkeypatch.setattr(georeferencer, "orthocorrection", mock_orthocorrection)
    dataset = reader.get_calibrated_dataset()

    original_lons = dataset["longitude"].copy()
    original_lats = dataset["latitude"].copy()

    dataset = georeferencer.orthocorrection(dataset)
    assert np.array_equal(dataset["longitude"], original_lons)
    assert np.array_equal(dataset["latitude"], original_lats)

    assert "tc_lons" in dataset
    assert "tc_lats" in dataset

    assert dataset["tc_lons"].shape == dataset["longitude"].shape
    assert dataset["tc_lats"].shape == dataset["latitude"].shape

    assert not np.array_equal(dataset["tc_lons"], dataset["longitude"])
    assert not np.array_equal(dataset["tc_lats"], dataset["latitude"])

def test_read_tle_file(pod_tle, tmp_path):
    tle_content = """1 33591U 09005A   24076.18425395  .00000218  00000+0  14176-3 0  9993
2 33591  99.0596 130.9575 0013809 190.5723 169.5160 14.12946284778493
1 33591U 09005A   24077.17564174  .00000205  00000+0  13478-3 0  9992
2 33591  99.0594 131.9606 0013864 187.8089 172.2868 14.12946798778634

"""
    tle_file = tmp_path / "test.tle"
    tle_file.write_text(tle_content)

    reader = LACPODReader(
        tle_dir=pod_tle.parent,
        tle_name=pod_tle.name,
        compute_lonlats_from_tles=True,
    )
    result = reader.read_tle_file(tle_file)

    expected = [
        "1 33591U 09005A   24076.18425395  .00000218  00000+0  14176-3 0  9993\n",
        "2 33591  99.0596 130.9575 0013809 190.5723 169.5160 14.12946284778493\n",
        "1 33591U 09005A   24077.17564174  .00000205  00000+0  13478-3 0  9992\n",
        "2 33591  99.0594 131.9606 0013864 187.8089 172.2868 14.12946798778634\n",
    ]
    assert result == expected


def test_missing_tle_file_raises_no_tle_data(pod_file_with_tbm_header, pod_tle):
    """A platform with no TLE file at all is missing orbit data, not a broken install.

    TLE_metopc.txt does not exist in the shipped bundle, so every Metop-C pass
    died with a bare FileNotFoundError that no caller was expecting, instead of
    the domain error the rest of the chain already handles.
    """
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name="TLE_%(satname)s_absent.txt",
                          compute_lonlats_from_tles=True)
    reader.read(pod_file_with_tbm_header)
    with pytest.raises(NoTLEData):
        reader.get_tle_lines()


def test_calibration_flag_does_not_void_geolocation(pod_file_with_tbm_header, pod_tle):
    """"Insufficient data for calibration" says nothing about where the pixel is.

    Two early NOAA-19 passes carry this bit on 96.8 % of their scan lines. It was
    folded into the mask applied to longitude and latitude, so almost the whole
    swath was voided even though the geolocation was sound - which then left
    omerc_bb without a usable nadir track.
    """
    from pygac.pod_reader import POD_QualityIndicator

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True)
    reader.read(pod_file_with_tbm_header)
    reader.scans["quality_indicators"] = POD_QualityIndicator.CALIBRATION

    lons, lats = reader.get_lonlat()
    assert np.isfinite(lons).any(), "calibration flag must not void the geolocation"

    # but a flag that does concern geolocation still must
    reader.lons = reader.lats = None
    reader._mask = None
    reader.scans["quality_indicators"] = POD_QualityIndicator.NO_EARTH_LOCATION
    lons, lats = reader.get_lonlat()
    assert not np.isfinite(lons).any(), "missing earth location must still void the geolocation"


class ReaderLackingClockOffsets(FakeGACReader):
    """A reader whose platform no clock drift table covers.

    ``pygac.clock_offsets_converter`` holds tables for noaa7, noaa9, noaa11,
    noaa12 and noaa14 only. For any other platform ``get_offsets`` raises a
    ``KeyError``, which is what this stands in for.
    """

    def compute_clock_offsets(self):
        """Refuse, the way the real reader does for an untabulated platform."""
        raise KeyError(self.spacecraft_name)

    def _compute_lonlat_from_tles(self, times):
        """Return coordinates of the right shape without needing an orbit."""
        return (np.zeros((len(times), 51)), np.zeros((len(times), 51)))


def test_a_platform_without_clock_offsets_is_still_geolocated():
    """Holding no clock drift table must cost the correction, not the whole pass."""
    reader = ReaderLackingClockOffsets(interpolate_coords=False,
                                       compute_lonlats_from_tles=True,
                                       adjust_clock_drift=True)
    lons, lats = reader.get_lonlat()
    assert lons.shape[0] == reader.along_track


def test_failed_pre_alignment_does_not_abandon_georeferencing(pod_file_with_tbm_header, pod_tle,
                                                               monkeypatch):
    """The coarse time fit is a pre-alignment; failing it must not discard the registration.

    It is only there to centre the swath before matching. When its optimiser
    reported non-convergence the exception propagated out of _georeference_data
    and the GCP stage never ran at all, so passes whose imagery was perfectly
    matchable kept uncorrected navigation.
    """
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    import pyorbital.geoloc_avhrr
    def refuse(*args, **kwargs):
        raise RuntimeError("Time offset estimation did not converge")
    monkeypatch.setattr(pyorbital.geoloc_avhrr, "estimate_time_offset", refuse)

    def mock_disp(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        return 0, (0, 0, 0), ([10000] * 60, [1000] * 60)
    use_a_stub_georeferencer(monkeypatch, mock_disp)

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True, reference_image="some_world_image.tif",
                          adjust_clock_drift=False)
    reader.read(pod_file_with_tbm_header)
    dataset = reader.get_calibrated_dataset()
    assert dataset.attrs["georeferenced"] is True
    assert dataset.attrs["median_gcp_distance"] == 1000


def test_the_navigation_travels_as_one_record(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """Everything the fit settled is gathered under a single attribute.

    Each hop between here and the delivered file forwards attributes by name, and
    each keeps its own list of the names it knows. A loose name therefore has to be
    added to every list along the way, and a renamed one is dropped in silence --
    which is how the fitted attitude came to be missing from every product. One name
    that never changes cannot go missing that way, and the mapping under it can grow
    without any hop being told.
    """
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    attitude = (0.001, -0.002, 0.003)   # radians, as pyorbital returns them

    def mock_disp(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        return 0, attitude, ([10000] * 60, [1000] * 60)
    use_a_stub_georeferencer(monkeypatch, mock_disp)

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True, reference_image="some_world_image.tif")
    reader.read(pod_file_with_tbm_header)
    dataset = reader.get_calibrated_dataset()

    roll, pitch, yaw = np.rad2deg(attitude)
    assert dataset.attrs["navigation"]["attitude_in_degrees"] == {
        "roll": pytest.approx(roll), "pitch": pytest.approx(pitch), "yaw": pytest.approx(yaw)}


def test_the_navigation_record_says_which_version_of_itself_it_is(pod_file_with_tbm_header,
                                                                 pod_tle, monkeypatch):
    """A consumer has to be able to tell what the record it is holding contains.

    The set of things recorded here has changed before and will change again, and a
    reader of an old product cannot ask the code that wrote it. The version travels
    inside the record rather than beside it, so that it cannot be separated from
    what it describes on the way to the file.
    """
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    def mock_disp(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        return 0, (0.001, 0.0, 0.003), ([10000] * 60, [1000] * 60)
    use_a_stub_georeferencer(monkeypatch, mock_disp)

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True, reference_image="some_world_image.tif")
    reader.read(pod_file_with_tbm_header)
    dataset = reader.get_calibrated_dataset()

    assert dataset.attrs["navigation"]["schema_version"] == NAVIGATION_METADATA_SCHEMA_VERSION


def test_the_navigation_record_carries_the_time_the_fit_settled(pod_file_with_tbm_header,
                                                                pod_tle, monkeypatch):
    """The time offset belongs with the attitude, under the one name that travels.

    It is half of the same answer: the fit divides an along-track displacement
    between a time and a pitch, and a consumer reading one without the other cannot
    tell which was fitted and which was held.
    """
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    fitted = 1.75

    def mock_disp(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        return fitted, (0.001, 0.0, 0.003), ([10000] * 60, [1000] * 60)
    use_a_stub_georeferencer(monkeypatch, mock_disp)

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True, reference_image="some_world_image.tif")
    reader.read(pod_file_with_tbm_header)
    dataset = reader.get_calibrated_dataset()

    assert dataset.attrs["navigation"]["time_offset_in_seconds"] == pytest.approx(fitted)


def test_estimated_attitude_is_published_in_degrees(pod_file_with_tbm_header, pod_tle,
                                                    monkeypatch):
    """The fitted attitude reaches the product in degrees.

    pyorbital returns the minimiser's own variables in radians, and pygac once
    stored them unconverted under a name claiming degrees, so anything screening
    on that field read a number 57 times too small. The name was corrected
    first; the unit is corrected here, because degrees is what the attitude is
    quoted in everywhere it is compared -- the OSPO notices, the published
    platform constants, and the analysis document.
    """
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    attitude = (0.001, -0.002, 0.003)   # radians, as pyorbital returns them
    def mock_disp(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        return 0, attitude, ([10000] * 60, [1000] * 60)
    use_a_stub_georeferencer(monkeypatch, mock_disp)

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True, reference_image="some_world_image.tif")
    reader.read(pod_file_with_tbm_header)
    dataset = reader.get_calibrated_dataset()
    assert dataset.attrs["estimated_attitude_in_degrees"] == pytest.approx(np.rad2deg(attitude))
    assert "estimated_attitude_in_radians" not in dataset.attrs
    assert dataset.attrs["navigation_metadata_schema_version"] >= 3


def test_the_navigation_record_says_what_the_fit_rested_on(pod_file_with_tbm_header, pod_tle,
                                                            monkeypatch):
    """Whether the navigation was corrected at all travels with what the fit settled.

    An attitude means nothing without knowing whether it was applied. The two belong
    to the same record, so that a reader cannot pick up one and miss the other.
    """
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    def mock_disp(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        return 0, (0, 0, 0), ([10000] * 60, [1000] * 60)
    use_a_stub_georeferencer(monkeypatch, mock_disp)

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True, reference_image="some_world_image.tif",
                          adjust_clock_drift=False)
    reader.read(pod_file_with_tbm_header)
    dataset = reader.get_calibrated_dataset()

    assert dataset.attrs["navigation"]["georeferenced"] is True


def navigated_with_too_few_gcps(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """Navigate a pass whose fit rests on one control point, and hand back the result.

    The displacement fit solves for four parameters (time, roll, pitch, yaw).
    With one surviving GCP it is wildly underdetermined: the optimiser drives
    the single residual to ~0 and reports a sub-millimetre median that looks
    like a perfect registration. The reader refuses it and warns.
    """
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    def mock_disp(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        # one control point, and a residual that looks flawless
        return 0, (0, 0, 0), ([10000], [0.002])
    use_a_stub_georeferencer(monkeypatch, mock_disp)

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True, reference_image="some_world_image.tif")
    reader.read(pod_file_with_tbm_header)
    with pytest.warns(RuntimeWarning):
        return reader.get_calibrated_dataset()


def test_georeferencing_rejects_too_few_gcps(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """A fit resting on a handful of control points must be rejected."""
    dataset = navigated_with_too_few_gcps(pod_file_with_tbm_header, pod_tle, monkeypatch)
    assert dataset.attrs["georeferenced"] is False


def test_a_reader_is_released_once_its_pass_is_done(pod_file_with_tbm_header, pod_tle):
    """Nothing may outlive the pass it belongs to.

    A cache on a bound method keys on the instance, so it keeps every reader that has
    ever run -- and the whole scanline array each one owns -- for the life of the
    process. Over a campaign that is measured in hundreds of thousands of passes, a
    worker that never lets go is a worker that eventually dies.
    """
    import gc
    import weakref

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True)
    reader.read(pod_file_with_tbm_header)
    reader.get_times()
    still_there = weakref.ref(reader)

    del reader
    gc.collect()

    assert still_there() is None


def test_a_reader_that_read_its_telemetry_is_released_too(pod_file_with_tbm_header, pod_tle):
    """The timestamps are not the only thing a bound cache holds on to."""
    import gc
    import weakref

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True)
    reader.read(pod_file_with_tbm_header)
    reader.get_telemetry()
    still_there = weakref.ref(reader)

    del reader
    gc.collect()

    assert still_there() is None


def test_the_corrupt_mask_can_be_asked_to_cover_more(pod_file_with_tbm_header, pod_tle):
    """A campaign may need to discard scanlines the library's default keeps.

    Pseudo-noise lines look perfectly well behaved along the scan -- one NOAA-14 pass
    agrees 0.999 between neighbouring pixels -- while being ruined down the track. The
    file flags them, and the default mask does not read that flag, so the matcher is
    handed every one of them.
    """
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True,
                          corrupt_flags="FATAL_FLAG|PSEUDO_NOISE")
    reader.read(pod_file_with_tbm_header)
    reader.scans["quality_indicators"] = POD_QualityIndicator.PSEUDO_NOISE

    assert reader.mask.all()


def test_a_misspelt_quality_flag_is_refused(pod_file_with_tbm_header, pod_tle):
    """A name that is not a quality flag must stop the pass, not quietly do nothing.

    The mask is named in configuration, so a typo would otherwise switch off part of it
    for every pass of a campaign without anything saying so.
    """
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True,
                          corrupt_flags="FATAL_FLAG|PSEUDONOISE")
    reader.read(pod_file_with_tbm_header)

    with pytest.raises(ReaderError, match="PSEUDONOISE"):
        reader.mask


def test_a_pass_that_could_not_be_calibrated_is_refused(pod_file_with_tbm_header, pod_tle):
    """A pass whose scanlines carry the calibration flag has no usable measurements left.

    Two NOAA-19 passes in the sample set it on 97% of their scanlines; the counts survive
    but the calibrated channels do not, leaving 3% of the imagery. Such a pass should be
    refused for what is wrong with it, rather than limping on to fail the point count.
    """
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True)
    reader.read(pod_file_with_tbm_header)
    reader.scans["quality_indicators"] = POD_QualityIndicator.CALIBRATION

    with pytest.raises(ReaderError, match="calibrated"):
        reader.get_calibrated_dataset()


def packed_into_words(counts):
    """Pack a (line, pixel, channel) array of counts the way the instrument does.

    Three ten-bit samples travel in each thirty-two bit word, taken in order across the
    channels of one pixel and on to the next, so a test that wants one channel to differ
    from another has to lay them out this way.
    """
    flat = counts.reshape(counts.shape[0], -1).astype(np.uint32)
    words = np.zeros((counts.shape[0], -(-flat.shape[1] // 3)), dtype=">u4")
    for offset, shift in ((0, 20), (1, 10), (2, 0)):
        taken = flat[:, offset::3]
        words[:, :taken.shape[1]] |= taken << shift
    return words


def test_a_pass_is_kept_when_any_one_channel_still_carries_an_image(pod_file_with_tbm_header,
                                                                    pod_tle, monkeypatch):
    """A dead visible channel does not condemn a pass whose thermal channel is sound.

    Matching falls back to channel 4 by night already, and there is nothing stopping it
    doing so by day, so a pass should only be refused when no channel holds an image.
    """
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True)
    reader.read(pod_file_with_tbm_header)
    lines, width = len(reader.scans), reader.scan_width
    counts = np.zeros((lines, width, 5), dtype=np.uint32)
    rng = np.random.default_rng(0)
    ramp = np.linspace(100, 900, width).astype(np.uint32)
    for channel in range(5):
        counts[:, :, channel] = rng.integers(0, 1024, size=(lines, width))
    counts[:, :, 3] = ramp                      # channel 4 holds a clean scene
    reader.scans["sensor_data"] = packed_into_words(counts)

    reader.get_calibrated_dataset()              # must not raise


def test_a_pass_carrying_no_image_is_refused(pod_file_with_tbm_header, pod_tle):
    """Noise is not a scene, however bright it is.

    Three NOAA-8 passes in the sample carry no coherent image at all: their neighbouring
    pixels agree no better than chance, where every pass that registers agrees above 0.89.
    Brightness cannot tell the two apart, since a snowfield is allowed to be bright.
    """
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True)
    reader.read(pod_file_with_tbm_header)
    noise = np.random.default_rng(0).integers(0, 1024, size=reader.scans["sensor_data"].shape)
    reader.scans["sensor_data"] = noise

    with pytest.raises(ReaderError, match="coherent"):
        reader.get_calibrated_dataset()


def test_georeferencing_accepts_a_dozen_control_points(pod_file_with_tbm_header, pod_tle,
                                                       monkeypatch):
    """Four parameters need two points for rank, so a dozen is a real constraint.

    The floor is there to catch the absurd one- and two-point "successes", not to insist
    on a rich harvest: a pass over ocean with one coastline in view can be well determined
    on a dozen well spread points, and passes that are not are refused on other grounds.
    """
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    def mock_disp(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds, count=12)
        return 0, (0, 0, 0), ([10000] * 12, [1000] * 12)
    use_a_stub_georeferencer(monkeypatch, mock_disp)

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True, reference_image="some_world_image.tif")
    reader.read(pod_file_with_tbm_header)
    dataset = reader.get_calibrated_dataset()

    assert dataset.attrs["georeferenced"] is True


def test_rejected_georeferencing_still_records_diagnostics(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """A rejected fit must still say how it was judged.

    The residual used to be written only after the quality gate passed, so a
    rejected pass carried no value at all and was indistinguishable downstream
    from one where registration was never attempted.
    """
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    def mock_disp(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds, count=8)
        return 0, (0, 0, 0), ([10000] * 8, [10000] * 8)   # rejected: too few points
    use_a_stub_georeferencer(monkeypatch, mock_disp)

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True, reference_image="some_world_image.tif")
    reader.read(pod_file_with_tbm_header)
    with pytest.warns(RuntimeWarning):
        dataset = reader.get_calibrated_dataset()
    assert dataset.attrs["georeferenced"] is False
    assert dataset.attrs["gcp_count"] == 8
    assert dataset.attrs["median_gcp_distance"] == 10000


def test_georeferencing_rejects_non_finite_residual(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """A NaN residual must be rejected, not accepted.

    `nan > 5000` is False, so a non-finite median used to sail through the
    quality gate. The fit that produced it had every attitude parameter pinned
    to its bound, and recomputing the geolocation from it voided the whole
    swath, losing the pass entirely.
    """
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    def mock_disp(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        return 0, (0, 0, 0), ([10000] * 60, [np.nan] * 60)
    use_a_stub_georeferencer(monkeypatch, mock_disp)

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True, reference_image="some_world_image.tif")
    reader.read(pod_file_with_tbm_header)
    with pytest.warns(RuntimeWarning):
        dataset = reader.get_calibrated_dataset()
    assert dataset.attrs["georeferenced"] is False


def test_georeferencing_rejects_attitude_on_its_bound(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """An attitude sitting on its optimiser bound is not a solution.

    The fit is bounded at +-0.5 rad, which is +-28.6 degrees; a real attitude
    error is a fraction of a degree. All three parameters resting on the bound
    means the minimiser had nothing to constrain it and walked into the corner
    of the box.
    """
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    def mock_disp(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        # plenty of points and a plausible residual, but a cornered attitude
        return 0, (0.5, -0.5, 0.5), ([10000] * 60, [1000] * 60)
    use_a_stub_georeferencer(monkeypatch, mock_disp)

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True, reference_image="some_world_image.tif")
    reader.read(pod_file_with_tbm_header)
    with pytest.warns(RuntimeWarning):
        dataset = reader.get_calibrated_dataset()
    assert dataset.attrs["georeferenced"] is False


def test_georeferencing_rejects_a_displacement_field_that_hangs_together_badly(
        pod_file_with_tbm_header, pod_tle, monkeypatch):
    """Matches that agree with no geometry are not a registration.

    Every parameter the fit solves for moves the swath smoothly, so a real
    displacement field is smooth too. One that is not was assembled from matches
    that found whatever the covariance surface happened to peak on, and fitting
    four parameters to it produces a number rather than a position.
    """
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    def mock_disp(calibrated_ds, *args, **rest):
        rng = np.random.default_rng(0)
        scattered = rng.uniform(-24, 24, size=60)
        for name, values in (("gcp_y", np.arange(60.0) * 30),
                             ("gcp_x", np.arange(60.0) * 30),
                             ("gcp_y_displacement", scattered),
                             ("gcp_x_displacement", rng.uniform(-24, 24, size=60))):
            calibrated_ds[name] = xr.DataArray(values, dims=["points"])
        return 0, (0, 0, 0), ([10000] * 60, [1000] * 60)
    use_a_stub_georeferencer(monkeypatch, mock_disp)

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True, reference_image="some_world_image.tif")
    reader.read(pod_file_with_tbm_header)
    with pytest.warns(RuntimeWarning):
        dataset = reader.get_calibrated_dataset()
    assert dataset.attrs["georeferenced"] is False


def test_a_large_residual_alone_does_not_reject_a_registration(pod_file_with_tbm_header, pod_tle,
                                                               monkeypatch):
    """How far the swath had to move is not evidence that the move was wrong.

    The residual is a distance between matched points and the reference after the
    fit has done its work. It is large when the navigation it started from was
    poor, which is the case the registration exists to repair. Judged on it, a
    pass that matched only near nadir passes while the same pass matched across
    the whole swath fails, because the pixel is six times larger at the edge.
    Whether the matches hang together is the question worth asking, and the
    coherence check asks it.
    """
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    def mock_disp(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        return 0, (0, 0, 0), ([20000] * 60, [10000] * 60)
    use_a_stub_georeferencer(monkeypatch, mock_disp)

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True, reference_image="some_world_image.tif")
    reader.read(pod_file_with_tbm_header)
    dataset = reader.get_calibrated_dataset()

    assert dataset.attrs["georeferenced"] is True
    assert dataset.attrs["median_gcp_distance"] == 10000


def test_metop_holds_its_swath_square_to_the_ground_track():
    """Metop turns as it flies, so its scan stays perpendicular to the ground track."""
    assert yaw_steers("metopa")


def test_an_unlisted_platform_scans_the_nominal_angle():
    """Appendix J steps the scan 0.05409868099 deg, putting sample 2048 at 55.37."""
    assert max_scan_angle_for("noaa14") == 55.37


def test_a_platform_measured_to_scan_differently_gets_its_own_angle(monkeypatch):
    """Metop-A, -B and -C were built together and scan 0.05 degrees apart."""
    monkeypatch.setitem(MAX_SCAN_ANGLES, "metopa", 55.245)
    assert max_scan_angle_for("metopa") == 55.245


def test_the_recorded_scan_angle_is_the_one_its_platform_scans(pod_file_with_tbm_header,
                                                               pod_tle, monkeypatch):
    """The attitude fit reads this attribute, so it has to name the swath it navigates."""
    reader = LACPODReader(interpolate_coords=True, tle_dir=pod_tle.parent, tle_name=pod_tle.name)
    reader.read(pod_file_with_tbm_header)
    monkeypatch.setitem(MAX_SCAN_ANGLES, reader.spacecraft_name, 55.245)

    assert reader.create_counts_dataset().attrs["max_scan_angle"] == 55.245


def test_the_dataset_carries_the_times_as_the_file_stated_them(pod_file_with_tbm_header, pod_tle):
    """A user of the product must be able to undo what pygac did to the timestamps.

    The scanline times are moved by the buffered-frame correction, by the clock
    drift adjustment and by the navigation fit, all in place. The product
    therefore carries what the file said as well as what pygac made of it.
    """
    reader = LACPODReader(interpolate_coords=True, tle_dir=pod_tle.parent, tle_name=pod_tle.name)
    reader.read(pod_file_with_tbm_header)

    dataset = reader.create_counts_dataset()

    numpy.testing.assert_array_equal(dataset["times_as_the_file_stated_them"],
                                     reader.times_as_the_file_stated_them)


def test_noaa16_records_the_angle_it_is_navigated_at(pod_file_with_tbm_header, pod_tle):
    """Its fit was handed a swath 0.12 degrees narrower than the one its navigation drew."""
    reader = LACPODReader(interpolate_coords=True, tle_dir=pod_tle.parent, tle_name=pod_tle.name)
    reader.read(pod_file_with_tbm_header)
    reader.spacecraft_name = "noaa16"

    recorded = reader.create_counts_dataset().attrs["max_scan_angle"]

    assert recorded == max_scan_angle_for("noaa16")


def test_a_platform_is_navigated_at_the_angle_it_scans(pod_file_with_tbm_header, pod_tle,
                                                      monkeypatch):
    """The swath the fit is given and the swath the navigation draws must be the same one."""
    def lonlats_scanning(angle):
        reader = LACPODReader(interpolate_coords=True, tle_dir=pod_tle.parent,
                              tle_name=pod_tle.name, compute_lonlats_from_tles=True)
        reader.read(pod_file_with_tbm_header)
        monkeypatch.setitem(MAX_SCAN_ANGLES, reader.spacecraft_name, angle)
        return reader.get_lonlat()[0]

    nominal = lonlats_scanning(55.37)
    narrower = lonlats_scanning(55.20)

    assert np.abs(narrower - nominal).max() > 0.01


def test_only_metop_is_credited_with_a_trusted_clock():
    """Every platform but Metop drifts at some point in its life.

    A POD platform drifts between resets, and a KLM platform's daily correction
    is not evidenced across the whole record. Crediting one of them with a
    trusted clock pins its time at zero and hands the residual to the pitch,
    which is the degeneracy this campaign keeps falling into. Metop is the one
    series whose clock is steered in flight throughout.
    """
    assert has_a_trusted_clock("metopa")
    assert has_a_trusted_clock("metopb")
    assert has_a_trusted_clock("metopc")
    assert not has_a_trusted_clock("noaa14")
    assert not has_a_trusted_clock("noaa19")
    assert not has_a_trusted_clock("noaa10")


def test_a_shift_of_scanlines_is_read_as_a_duration():
    """The georeferencer counts scanlines; only pygac knows how fast they arrive.

    Turning a count of lines into a duration needs the scan rate, which is a fact
    about the instrument. Keeping that conversion here leaves the georeferencer
    free of any assumption about what it is looking at, so it can be pointed at
    GAC or at another instrument without change.
    """
    six_per_second = np.timedelta64(166667, "us")
    times = np.datetime64("1997-11-09T19:00:00") + np.arange(100) * six_per_second

    assert seconds_from_scanlines(163, times) == pytest.approx(27.17, abs=0.01)


def test_the_early_klm_frame_buffer_lag_is_dated_and_per_platform():
    """NOAA-15 and NOAA-16 stamped buffered scanlines with the time of readout.

    The onboard software held eleven frames where five were assumed, so the data
    handed out was six frames -- 900 ms at 150 ms a frame -- older than the time
    written beside it. NOAA corrected it on 2001-08-07 at 23:59Z. This is a fault
    in what the platform reported, attested outside this record, so it is dated
    rather than fitted.

    Measured on the six early NOAA-15 and NOAA-16 passes we hold: taking the lag
    out moves the mean fitted time offset from 1.367 s to 0.745 s, while adding
    it instead pushes them to 2.556 s.
    """
    before = np.datetime64("2001-06-01T00:00:00")
    after = np.datetime64("2001-09-01T00:00:00")

    assert frame_buffer_lag_for("noaa15", before) == 0.9
    assert frame_buffer_lag_for("noaa16", before) == 0.9
    assert frame_buffer_lag_for("noaa15", after) == 0.0
    assert frame_buffer_lag_for("noaa16", after) == 0.0
    assert frame_buffer_lag_for("noaa17", before) == 0.0
    assert frame_buffer_lag_for("noaa14", before) == 0.0


def test_a_drifting_platform_has_its_pitch_held_while_its_time_is_fitted(
        pod_file_with_tbm_header, pod_tle, monkeypatch):
    """The policy has to reach the fit, not merely be stated beside it.

    noaa14's clock is not believed, so its time is fitted and its pitch is the one
    held at zero. A fit told to solve for both would hand back a pitch invented out
    of the same displacement the time already accounts for.
    """
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    asked = {}

    def measures_no_drift(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        return np.zeros((60, 2)), np.zeros((60, 2)), 0.0

    def fits(gcps, ref_lons, ref_lats, start_time, tle, max_scan_angle, **rest):
        asked.update(rest)
        return 0.0, (0, 0, 0), ([10000] * 60, [1000] * 60)

    from georeferencer import georeferencer
    from pyorbital import geoloc_avhrr
    monkeypatch.setattr(georeferencer, "measure_swath_displacement", measures_no_drift)
    monkeypatch.setattr(geoloc_avhrr, "fit_navigation", answering_as_a_navigation_fit(fits))

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True, reference_image="some_world_image.tif")
    reader.read(pod_file_with_tbm_header)
    reader.spacecraft_name = "noaa14"
    reader.get_calibrated_dataset()

    assert asked["solve_for_time"] is True
    assert asked["solve_for_pitch"] is False


def test_exactly_one_of_the_time_and_the_pitch_is_ever_fitted():
    """The two say the same thing about a pass, so fitting both leaves both meaningless.

    A shift along the track can be written as a clock error or as a pitch, and they
    exchange at about 2.25 seconds per degree. Only the curvature the shift leaves
    across the swath separates them, which a single pass rarely constrains, so a fit
    given both will divide the shift between them arbitrarily and report two numbers
    where the data supports one. Whichever is not fitted is held at zero.
    """
    for platform in ("metopa", "metopb", "metopc", "noaa19", "noaa16", "noaa14", "noaa10", "noaa7"):
        assert should_fit_the_clock(platform) != should_fit_the_pitch(platform)


def test_whether_to_fit_the_clock_is_settled_by_the_platform_alone():
    """The choice is a property of the platform and of nothing about the pass.

    Deciding per pass -- on how far along its track a pass sits, or on whether a
    clock table reaches it -- lets the same platform report a time offset on one
    pass and a pitch on the next, from data that cannot tell the two apart. Every
    platform but Metop therefore has its time fitted, always.
    """
    assert should_fit_the_clock("noaa10")
    assert should_fit_the_clock("noaa14")
    assert should_fit_the_clock("noaa19")
    assert not should_fit_the_clock("metopa")
    assert not should_fit_the_clock("metopb")


def test_a_pass_past_its_platform_table_is_fitted_end_to_end(pod_file_with_tbm_header,
                                                            pod_tle, monkeypatch):
    """The fixture is a NOAA-14 pass from November 2000 and the table stops that August,
    so nothing measured reaches it however still the swath sits."""
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    fitted = 2.5

    def measures_no_drift(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        return np.zeros((60, 2)), np.zeros((60, 2)), 0.0

    def fits(gcps, ref_lons, ref_lats, start_time, tle, max_scan_angle, solve_for_time=True, **rest):
        return (fitted if solve_for_time else 0.0), (0, 0, 0), ([10000] * 60, [1000] * 60)

    from georeferencer import georeferencer
    from pyorbital import geoloc_avhrr
    monkeypatch.setattr(georeferencer, "measure_swath_displacement", measures_no_drift)
    monkeypatch.setattr(geoloc_avhrr, "fit_navigation", answering_as_a_navigation_fit(fits))
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True, reference_image="some_world_image.tif")
    reader.read(pod_file_with_tbm_header)
    reader.spacecraft_name = "noaa14"

    dataset = reader.get_calibrated_dataset()

    assert dataset.attrs["estimated_time_offset_in_seconds"] == pytest.approx(fitted)


def test_a_pass_records_whether_a_clock_measurement_reached_it(pod_file_with_tbm_header,
                                                              pod_tle, monkeypatch):
    """Whether the time was fitted or taken from a table is provenance the campaign has
    to be able to read back off the product."""
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    def measures_no_drift(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        return np.zeros((60, 2)), np.zeros((60, 2)), 0.0

    def fits(gcps, ref_lons, ref_lats, start_time, tle, max_scan_angle, **rest):
        return 0.0, (0, 0, 0), ([10000] * 60, [1000] * 60)

    from georeferencer import georeferencer
    from pyorbital import geoloc_avhrr
    monkeypatch.setattr(georeferencer, "measure_swath_displacement", measures_no_drift)
    monkeypatch.setattr(geoloc_avhrr, "fit_navigation", answering_as_a_navigation_fit(fits))
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True, reference_image="some_world_image.tif")
    reader.read(pod_file_with_tbm_header)
    reader.spacecraft_name = "noaa14"

    dataset = reader.get_calibrated_dataset()

    assert dataset.attrs["clock_table_covers_the_pass"] is False


def test_noaa16_scans_the_narrower_angle_it_was_always_given():
    """Measured at 55.22, against 55.34 for the rest; the historical value stands until
    a wider sample than two consecutive days is in."""
    assert max_scan_angle_for("noaa16") == 55.25


def test_a_scan_angle_named_in_the_configuration_is_used(tmp_path):
    """An operator retunes a platform without editing pygac, as an FDR runner must allow."""
    from pygac.configuration import read_config_file, reset_config
    settings = tmp_path / "pygac.cfg"
    settings.write_text("[scan_angles]\nnoaa19 = 55.301\n")
    read_config_file(str(settings))
    try:
        assert max_scan_angle_for("noaa19") == 55.301
    finally:
        reset_config()


def test_a_scan_angle_is_read_from_the_named_configuration_file(tmp_path, monkeypatch):
    """A caller that drives the reader directly still gets the configured angle.

    Naming the file in ``PYGAC_CONFIG_FILE`` is how pygac has always been pointed at
    its settings, and it is the only handle a runner embedding the reader has. Only
    pygac's own command-line runner ever opened that file, so everything else -- the
    FDR chain among them -- was quietly navigated at the built-in angle while a
    configured one sat unread beside it.
    """
    from pygac.configuration import reset_config
    settings = tmp_path / "pygac.cfg"
    settings.write_text("[scan_angles]\nnoaa19 = 55.301\n")
    monkeypatch.setenv("PYGAC_CONFIG_FILE", str(settings))
    reset_config()
    try:
        assert max_scan_angle_for("noaa19") == 55.301
    finally:
        reset_config()


def test_the_navigation_record_says_how_many_points_it_rested_on(pod_file_with_tbm_header,
                                                                 pod_tle, monkeypatch):
    """An attitude means nothing without the control points behind it.

    The same angles read very differently resting on a dozen points and on six
    hundred, so the count belongs with them rather than beside them. It is the fact
    that went missing from delivered products: recorded on the dataset, dropped at
    the next hop, and absent from every file the campaign wrote.
    """
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    def mock_disp(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        return 0, (0, 0, 0), ([10000] * 60, [1000] * 60)
    use_a_stub_georeferencer(monkeypatch, mock_disp)

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True, reference_image="some_world_image.tif",
                          adjust_clock_drift=False)
    reader.read(pod_file_with_tbm_header)
    dataset = reader.get_calibrated_dataset()

    assert dataset.attrs["navigation"]["gcp_count"] == 60


def test_the_counts_are_carried_at_the_width_they_need(pod_file_with_tbm_header, pod_tle):
    """Ten-bit counts do not need sixty-four bits to hold them.

    Everything downstream inherits the width the counts are read at: the calibrated
    channels, the copy kept for the uncertainties, and every array derived from them.
    A pass carries two of these cubes, so the choice is worth hundreds of megabytes
    on a record of nine hundred thousand passes.
    """
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=False, adjust_clock_drift=False)
    reader.read(pod_file_with_tbm_header)

    assert reader.get_counts().dtype == np.float32


def test_the_navigation_record_says_how_old_its_elements_were(pod_file_with_tbm_header,
                                                               pod_tle, monkeypatch):
    """How far the orbit was propagated belongs with what was navigated from it.

    An element drifts about a kilometre a day, so a pass navigated from a week-old
    element starts several kilometres out along its own track. The fit absorbs most of
    that as a time offset, which is exactly why the age must travel with the product:
    without it there is no way to tell a clock error from an orbit propagated too far,
    and no way to weigh a fitted offset when it is read back as a clock measurement.
    """
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)

    def mock_disp(calibrated_ds, *args, **rest):
        record_a_coherent_field(calibrated_ds)
        return 0, (0, 0, 0), ([10000] * 60, [1000] * 60)
    use_a_stub_georeferencer(monkeypatch, mock_disp)

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True, reference_image="some_world_image.tif",
                          adjust_clock_drift=False)
    reader.read(pod_file_with_tbm_header)
    dataset = reader.get_calibrated_dataset()

    # the fixture's first scanline is 2000-11-17T03:18:44.8 and the nearer element's
    # epoch is 01:07:52.4, which is where the age comes from -- not the start time the
    # data set name states, which the reader does not use
    assert dataset.attrs["navigation"]["element_age_in_days"] == pytest.approx(0.0909, abs=0.0005)


def test_a_clock_named_in_the_configuration_is_believed_or_not(tmp_path):
    """A platform shown to hold its time is promoted without editing pygac.

    Whether a clock can be navigated on is a finding about a platform, not a fact
    about this code, so an FDR runner must be able to state it per platform.
    """
    from pygac.configuration import read_config_file, reset_config
    settings = tmp_path / "pygac.cfg"
    settings.write_text("[trusted_clocks]\nnoaa19 = yes\nmetopb = no\n")
    read_config_file(str(settings))
    try:
        assert has_a_trusted_clock("noaa19")
        assert not has_a_trusted_clock("metopb")
    finally:
        reset_config()


def test_the_poes_platforms_fly_without_turning():
    """NOAA POES holds a fixed attitude, so its scan follows the inertial track."""
    assert not yaw_steers("noaa19")


def test_a_steered_platform_is_navigated_differently(pod_file_with_tbm_header, pod_tle):
    """The platform decides whether the geolocation turns with the ground track."""
    def lonlats_flying_as(platform):
        reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                              compute_lonlats_from_tles=True)
        reader.read(pod_file_with_tbm_header)
        reader.spacecraft_name = platform
        return reader.get_lonlat()

    straight, _ = lonlats_flying_as("noaa19")
    turned, _ = lonlats_flying_as("metopa")

    assert np.abs(turned - straight).max() > 0.1


def test_the_geolocation_names_its_rotation_order(pod_file_with_tbm_header, pod_tle):
    """The geolocation chooses an order rather than falling back to the legacy one.

    pyorbital says so itself: asked to apply an attitude without being told which order
    to apply it in, it warns that it is using the order it shipped with, which its own
    validation shows is the less accurate one. It stays quiet when the platform is not
    pitched, because the two orders are then identical, so the platform has to be
    pitched for the question to arise at all.
    """
    import warnings as warnings_module

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True)
    reader.read(pod_file_with_tbm_header)
    reader._rpy = np.array([0.0, 0.01, 0.0])

    with warnings_module.catch_warnings(record=True) as caught:
        warnings_module.simplefilter("always")
        reader.get_lonlat()

    assert not [one for one in caught if "rotation order" in str(one.message)]


def test_the_geolocation_names_its_nadir_convention(pod_file_with_tbm_header, pod_tle):
    """The geolocation chooses a convention rather than falling back to the legacy one."""
    import warnings as warnings_module

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name,
                          compute_lonlats_from_tles=True)
    reader.read(pod_file_with_tbm_header)

    with warnings_module.catch_warnings(record=True) as caught:
        warnings_module.simplefilter("always")
        reader.get_lonlat()

    assert not [one for one in caught if "nadir convention" in str(one.message)]


def test_a_refused_georeferencing_says_so_in_the_navigation_metadata(
        pod_file_with_tbm_header, pod_tle, monkeypatch):
    """A pass that failed to georeference must say so where the product records it.

    The flag beside the rest of the navigation metadata is the one a downstream
    reader consults. Leaving it behind while the top-level flag turns False hands
    out uncorrected geolocation labelled as registered.
    """
    dataset = navigated_with_too_few_gcps(pod_file_with_tbm_header, pod_tle, monkeypatch)
    assert dataset.attrs["navigation"]["georeferenced"] is False


def test_a_refused_georeferencing_names_the_pass_in_the_log(
        pod_file_with_tbm_header, pod_tle, monkeypatch, caplog):
    """A campaign log must say which pass lost its registration, not only that one did.

    A run processes hundreds of thousands of passes. A traceback naming no file
    cannot be attributed to one, so a pass that quietly fell back to uncorrected
    geolocation reads as noise rather than as work to redo.
    """
    with caplog.at_level(logging.ERROR, logger="pygac.reader"):
        navigated_with_too_few_gcps(pod_file_with_tbm_header, pod_tle, monkeypatch)

    assert os.path.basename(pod_file_with_tbm_header) in caplog.text


def test_a_comment_in_a_tle_file_is_not_taken_for_an_element(tmp_path):
    """A file may say which of its elements were made rather than observed.

    The lines are paired by position, so anything that is not part of an element
    shifts every element after it -- the first line of one read together with the
    second line of the next -- and nothing reports it.
    """
    marked = tmp_path / "noaa14.tle"
    marked.write_text(
        "# elements below are synthetic, fitted to NOAA TBUS bulletins\n"
        "1 23455U 94089A   00322.04713399  .00000318  00000-0  19705-3 0  5298\n"
        "2 23455  99.1591 303.5706 0010037  25.7760 334.3905 14.12496755303183\n"
        "# and this one came from the archive\n"
        "1 23455U 94089A   00322.96799836  .00000229  00000-0  14918-3 0  5303\n"
        "2 23455  99.1590 304.5117 0009979  23.1101 337.0518 14.12496633303313\n")

    reader = LACPODReader(tle_dir=str(tmp_path), tle_name="noaa14.tle")
    kept = reader.read_tle_file(str(marked))

    assert [line[0] for line in kept] == ["1", "2", "1", "2"]


def a_georeferenced_pass(pod_file_with_tbm_header, pod_tle, monkeypatch, platform=None, **reader_kwargs):
    """Read and calibrate the fixture pass with georeferencing on, whatever stands in for it."""
    def skip_thermal(channels, *args, **kwargs):
        return channels, []
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal_channels", skip_thermal)
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name, compute_lonlats_from_tles=True,
                          reference_image="some_world_image.tif", **reader_kwargs)
    reader.read(pod_file_with_tbm_header)
    if platform:
        reader.spacecraft_name = platform
    return reader


def test_every_control_point_keeps_what_the_fit_left_of_its_miss(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """A median says what a typical point does; how well each one was matched needs each one.

    The leftovers are kept in swath pixels, next to the displacements they are what is
    left of, so structure the fit could not absorb -- a drift along the track, a bend
    across it -- can be seen and fitted later.
    """
    misses = np.column_stack([np.linspace(-0.4, 0.4, 60), np.linspace(0.3, -0.3, 60)])
    misses[5] = np.nan
    navigate_with(monkeypatch, misses=misses)
    reader = a_georeferenced_pass(pod_file_with_tbm_header, pod_tle, monkeypatch)

    dataset = reader.get_calibrated_dataset()

    np.testing.assert_array_equal(dataset["gcp_y_residual"].values, misses[:, 0])
    np.testing.assert_array_equal(dataset["gcp_x_residual"].values, misses[:, 1])
    assert dataset["gcp_y_residual"].dims == ("points",)
    assert dataset["gcp_x_residual"].attrs["units"] == "swath pixels"


def test_a_refused_pass_keeps_its_leftovers_too(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """A pass turned away still has to show why, point by point."""
    navigate_with(monkeypatch, count=8, misses=np.full((8, 2), 0.25))
    reader = a_georeferenced_pass(pod_file_with_tbm_header, pod_tle, monkeypatch)

    with pytest.warns(RuntimeWarning):
        dataset = reader.get_calibrated_dataset()

    assert dataset.attrs["georeferenced"] is False
    np.testing.assert_array_equal(dataset["gcp_x_residual"].values, np.full(8, 0.25))


def test_the_spread_of_the_leftovers_is_summarised_along_and_across(pod_file_with_tbm_header, pod_tle,
                                                                    monkeypatch):
    """The spread is quoted per axis and with its tail, since a median hides both.

    Lines and samples are kept apart because they answer to different parameters:
    along the track to time and pitch, across it to roll and yaw. The summary joins
    the navigation record beside what the fit changed, rather than replacing it.
    """
    misses = np.column_stack([np.full(60, 0.3), np.tile([0.4, -0.4], 30)])
    misses[[0, 1]] = [[3.0, 4.0], [np.nan, np.nan]]
    navigate_with(monkeypatch, misses=misses)
    reader = a_georeferenced_pass(pod_file_with_tbm_header, pod_tle, monkeypatch)

    record = reader.get_calibrated_dataset().attrs["navigation"]

    assert record["residual_rms_in_lines"] == pytest.approx(np.sqrt((9.0 + 0.09 * 58) / 59))
    assert record["residual_rms_in_samples"] == pytest.approx(np.sqrt((16.0 + 0.16 * 58) / 59))
    assert record["residual_p95_in_pixels"] == pytest.approx(
        np.percentile(np.hypot(*np.vstack([[3.0, 4.0], misses[2:]]).T), 95))
    assert record["georeferenced"] is True
    assert "attitude_in_degrees" in record


def test_a_refused_pass_keeps_the_spread_in_its_record(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """A pass turned away is told apart from one never tried by what its record still says."""
    navigate_with(monkeypatch, count=8, misses=np.full((8, 2), 0.25))
    reader = a_georeferenced_pass(pod_file_with_tbm_header, pod_tle, monkeypatch)

    with pytest.warns(RuntimeWarning):
        record = reader.get_calibrated_dataset().attrs["navigation"]

    assert record["georeferenced"] is False
    assert record["residual_rms_in_samples"] == pytest.approx(0.25)


def test_the_record_says_how_far_its_navigation_can_be_trusted(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """The fitted time and attitude come with their uncertainty, in the units the record quotes.

    The record publishes the attitude in degrees, so its uncertainty and covariance are
    in degrees too. One of time and pitch is always held, and a held parameter was
    asserted rather than measured: it carries no uncertainty.
    """
    spread = np.array([0.1, 2e-4, 0.0, 4e-4])   # seconds and radians; the pitch was held
    correlation = np.eye(4)
    correlation[0, 1] = correlation[1, 0] = 0.3
    navigate_with(monkeypatch, covariance=correlation * np.outer(spread, spread))
    reader = a_georeferenced_pass(pod_file_with_tbm_header, pod_tle, monkeypatch)

    record = reader.get_calibrated_dataset().attrs["navigation"]

    in_degrees = spread * [1.0, *[180 / np.pi] * 3]
    assert record["time_offset_uncertainty_in_seconds"] == pytest.approx(0.1)
    assert record["attitude_uncertainty_in_degrees"] == {
        "roll": pytest.approx(in_degrees[1]), "pitch": pytest.approx(in_degrees[2]),
        "yaw": pytest.approx(in_degrees[3])}
    assert "time_pitch_correlation" not in record
    # a netCDF attribute holds one dimension, so the matrix is kept row by row, and says in what order
    assert record["covariance_order"] == "time_offset, roll, pitch, yaw"
    flat = record["covariance_in_seconds_and_degrees"]
    assert np.ndim(flat) == 1
    np.testing.assert_allclose(np.reshape(flat, (4, 4)), correlation * np.outer(in_degrees, in_degrees))


def test_the_record_says_what_its_navigation_was_fitted_with(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """Refitting a pass later is only comparable if every choice of the first fit can be reproduced.

    The coarse shift the matching found is kept too: the clock is decided by the platform
    alone, and how far along its track a pass sat is what would show that decision wrong.
    """
    from importlib.metadata import version

    navigate_with(monkeypatch, settings={"loss": "soft_l1", "solve_for_time": True}, along_track_lines=12)
    reader = a_georeferenced_pass(pod_file_with_tbm_header, pod_tle, monkeypatch)
    seconds_per_line = np.diff(reader.get_times()).mean() / np.timedelta64(1, "s")

    record = reader.get_calibrated_dataset().attrs["navigation"]

    assert record["fit_settings"] == {"loss": "soft_l1", "solve_for_time": True}
    assert record["pyorbital_version"] == version("pyorbital")
    assert record["georeferencer_version"] == version("georeferencer")
    assert record["coarse_along_track_shift_in_seconds"] == pytest.approx(12 * seconds_per_line)


def test_the_two_time_offsets_are_kept_apart(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """The fit starts from the pre-aligned times, so reproducing it needs both parts, not their sum."""
    navigate_with(monkeypatch, time_offset=0.25)
    reader = a_georeferenced_pass(pod_file_with_tbm_header, pod_tle, monkeypatch, adjust_clock_drift=False)
    monkeypatch.setattr(reader, "_correct_time_offset", lambda calibrated_ds: 1.5)

    record = reader.get_calibrated_dataset().attrs["navigation"]

    assert record["pre_alignment_time_offset_in_seconds"] == 1.5
    assert record["fitted_time_offset_in_seconds"] == 0.25
    assert record["time_offset_in_seconds"] == 1.75


def test_products_carrying_the_fit_quality_are_told_apart(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """Products written before the leftovers and uncertainties were kept must be distinguishable."""
    navigate_with(monkeypatch)
    reader = a_georeferenced_pass(pod_file_with_tbm_header, pod_tle, monkeypatch)

    dataset = reader.get_calibrated_dataset()

    assert dataset.attrs["navigation"]["schema_version"] == 7
    assert dataset.attrs["navigation_metadata_schema_version"] == 7


def test_each_pass_says_in_the_log_how_well_it_was_navigated(pod_file_with_tbm_header, pod_tle, monkeypatch, caplog):
    """A run log shows, pass by pass and without asking for debug output, how the fit went.

    The figures are the ones the product records: how many points the fit stood on,
    what it changed, and the spread of what it left in lines and samples. They are
    said before the fit is judged, so a refused pass shows them beside its reason.
    """
    import logging

    misses = np.column_stack([np.full(60, 0.3), np.full(60, -0.4)])
    navigate_with(monkeypatch, time_offset=0.44, misses=misses)
    reader = a_georeferenced_pass(pod_file_with_tbm_header, pod_tle, monkeypatch)

    with caplog.at_level(logging.INFO, logger="pygac.reader"):
        reader.get_calibrated_dataset()

    said = [record.getMessage() for record in caplog.records if record.levelno == logging.INFO]
    summary = [message for message in said if "navigation fit" in message]
    assert len(summary) == 1, said
    assert "60 control points" in summary[0]
    assert "time offset 0.44 s" in summary[0]
    assert "rms 0.30 lines, 0.40 samples" in summary[0]


def test_the_record_says_how_far_the_trustworthy_points_ended_up(pod_file_with_tbm_header, pod_tle, monkeypatch,
                                                                  caplog):
    """How good the corrected navigation is, apart from the matcher's noise and its mistakes.

    Mismatches are set aside and counted; the spread of the rest is quoted; and the
    part of it that varies smoothly across the swath is quoted on its own, as the
    navigation the fit did not remove. The pass's log line says the same.
    """
    import logging

    rng = np.random.default_rng(0)
    across = np.linspace(-1.0, 1.0, 300)
    misses = np.column_stack([np.zeros(300), 0.8 * across]) + rng.normal(0.0, 0.3, (300, 2))
    misses[[3, 50, 120]] += [[10.0, -14.0]] * 3
    navigate_with(monkeypatch, count=300, misses=misses)
    reader = a_georeferenced_pass(pod_file_with_tbm_header, pod_tle, monkeypatch)

    with caplog.at_level(logging.INFO, logger="pygac.reader"):
        record = reader.get_calibrated_dataset().attrs["navigation"]

    assert 3 <= record["residual_mismatch_count"] <= 4
    np.testing.assert_allclose(record["residual_trusted_rms_in_lines"], 0.3, rtol=0.2)
    assert record["residual_systematic_rms_in_samples"] == pytest.approx(0.8 * np.sqrt(np.mean(across ** 2)), rel=0.15)
    assert record["residual_systematic_rms_in_lines"] < 0.1
    assert "mismatches" in caplog.text and "smooth part" in caplog.text


def test_the_fit_is_told_the_terrain_height_of_every_control_point(pod_file_with_tbm_header, pod_tle, monkeypatch,
                                                                    tmp_path):
    """Given a DEM, the fit allows for how high each reference stands.

    A mountain is seen away from nadir of where it stands, and a fit not told so takes
    that for navigation error and tilts the whole pass to meet it. The heights are kept
    per point, next to the misses they explain.
    """
    import rasterio
    from rasterio.transform import from_origin

    dem = tmp_path / "dem.tif"
    heights = np.arange(10.0)[None, :].repeat(10, axis=0) * 500.0     # rising eastwards, 500 m a column
    with rasterio.open(dem, "w", driver="GTiff", height=10, width=10, count=1, dtype="float32", crs="EPSG:4326",
                       transform=from_origin(0.0, 10.0, 1.0, 1.0)) as out:
        out.write(heights.astype("float32"), 1)
    lonlats = np.column_stack([np.linspace(0.5, 9.5, 60), np.full(60, 5.0)])
    asked = []
    navigate_with(monkeypatch, gcp_lonlats=lonlats, fits=asked)
    from georeferencer import georeferencer
    monkeypatch.setattr(georeferencer, "orthocorrection", lambda calibrated_ds, *args: calibrated_ds)
    reader = a_georeferenced_pass(pod_file_with_tbm_header, pod_tle, monkeypatch, dem=str(dem))

    dataset = reader.get_calibrated_dataset()

    expected = np.floor(lonlats[:, 0]) * 500.0
    np.testing.assert_allclose(asked[0]["heights"], expected)
    np.testing.assert_allclose(dataset["gcp_height"].values, expected)
    assert dataset["gcp_height"].attrs["units"] == "m"


def test_without_a_dem_the_references_are_taken_at_sea_level(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """No DEM, no heights: the fit compares on the ellipsoid, as it always has."""
    asked = []
    navigate_with(monkeypatch, fits=asked)
    reader = a_georeferenced_pass(pod_file_with_tbm_header, pod_tle, monkeypatch)

    dataset = reader.get_calibrated_dataset()

    assert asked[0].get("heights") is None
    assert "gcp_height" not in dataset


def test_the_record_carries_how_the_pass_was_prepared_and_judged(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """The record is what reaches the product whole, so what the pass was judged on belongs in it.

    Whether the swath was pre-aligned, whether a clock measurement reached it, and how
    much of its displacement field no smooth geometry explains were left beside the
    record, and the reader forwarding the record dropped them on the way to the file.
    """
    navigate_with(monkeypatch)
    reader = a_georeferenced_pass(pod_file_with_tbm_header, pod_tle, monkeypatch, adjust_clock_drift=False)
    monkeypatch.setattr(reader, "_correct_time_offset", lambda calibrated_ds: 1.5)

    record = reader.get_calibrated_dataset().attrs["navigation"]

    assert record["pre_alignment_applied"] is True
    assert record["clock_table_covers_the_pass"] in (True, False)
    assert record["unexplained_displacement_in_pixels"] == pytest.approx(0.0, abs=1e-6)


def test_a_pass_refused_for_its_matches_keeps_what_they_were_judged_on(pod_file_with_tbm_header, pod_tle,
                                                                      monkeypatch):
    """A pass turned away for matches no geometry explains has to show how far they strayed."""
    def scattered(calibrated_ds, count):
        rng = np.random.default_rng(0)
        for name, values in (("gcp_y", np.arange(float(count)) * 30), ("gcp_x", np.arange(float(count)) * 30),
                             ("gcp_y_displacement", rng.uniform(-24, 24, count)),
                             ("gcp_x_displacement", rng.uniform(-24, 24, count))):
            calibrated_ds[name] = xr.DataArray(values, dims=["points"])
    navigate_with(monkeypatch, field=scattered)
    reader = a_georeferenced_pass(pod_file_with_tbm_header, pod_tle, monkeypatch)

    with pytest.warns(RuntimeWarning):
        record = reader.get_calibrated_dataset().attrs["navigation"]

    assert record["georeferenced"] is False
    assert record["unexplained_displacement_in_pixels"] > 8


def a_swath_crossing(lat, lon, heading, width, lines=120):
    """Return lon/lats laid out like a swath: lines along a track, samples spread across it."""
    from pyproj import Geod

    geod = Geod(ellps="WGS84")
    track_lons, track_lats, _ = geod.fwd(np.full(lines, lon), np.full(lines, lat),
                                         np.full(lines, heading), np.arange(lines) * 1100.0)
    across = np.tan(np.deg2rad(np.linspace(-55.37, 55.37, width))) * 830e3
    lons, lats, _ = geod.fwd(np.repeat(track_lons, width).reshape(lines, width),
                             np.repeat(track_lats, width).reshape(lines, width),
                             np.full((lines, width), heading + 90.0), np.tile(across, (lines, 1)))
    return lons, lats


@pytest.mark.parametrize("reader_class", ["GACPODReader", "GACKLMReader", "LACPODReader"])
@pytest.mark.parametrize(("lat", "lon", "heading"), [(-30.0, 20.0, -12.0), (70.0, 170.0, -5.0), (85.0, 0.0, 0.0)])
def test_the_tie_points_are_interpolated_along_the_satellite_spline(reader_class, lat, lon, heading):
    """Each scan line's pixels lie on geotiepoints' cubic spline through its tie points.

    The spline runs through earth-centred coordinates, so neither the dateline nor the
    pole bends it; past the first and last tie points it carries on as the same cubic.
    However it is computed, the pixels must land where that spline puts them.
    """
    import geotiepoints as gtp
    from pyproj import Geod

    import pygac.gac_klm
    import pygac.gac_pod
    import pygac.lac_pod

    modules = {"GACPODReader": pygac.gac_pod, "GACKLMReader": pygac.gac_klm, "LACPODReader": pygac.lac_pod}
    reader = getattr(modules[reader_class], reader_class)()
    lons, lats = a_swath_crossing(lat, lon, heading, reader.scan_width)
    ties = np.round(reader.lonlat_sample_points).astype(int)
    tie_lons, tie_lats = lons[:, ties], lats[:, ties]
    rows = np.arange(len(tie_lats))
    expected = gtp.SatelliteInterpolator((tie_lons, tie_lats), (rows, reader.lonlat_sample_points),
                                         (rows, np.arange(reader.scan_width)), 1, 3).interpolate()

    found = reader.lonlat_interpolator(tie_lons, tie_lats)

    _, _, miss = Geod(ellps="WGS84").inv(found[0].ravel(), found[1].ravel(), expected[0].ravel(), expected[1].ravel())
    assert miss.max() < 1e-3


def test_the_angles_are_computed_block_by_block_on_the_configured_dask_scheduler(pod_file_with_tbm_header,
                                                                                 pod_tle, monkeypatch):
    """The angles of a few scan lines at a time are dask tasks, run on the chain's scheduler.

    Every angle is worked out pixel by pixel, from the pixel's time and lon/lat alone,
    so taken a block of scan lines at a time, side by side, they come out as computed
    over the whole pass at once.
    """
    import dask
    from pyorbital import astronomy
    from pyorbital.orbital import Orbital

    import pygac.reader as reader_module
    from pygac.utils import centered_modulus, get_absolute_azimuth_angle_diff

    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name, compute_lonlats_from_tles=True)
    reader.read(pod_file_with_tbm_header)
    monkeypatch.setattr(reader_module, "ANGLE_LINES_PER_BLOCK", 1)
    computes = []

    def recording_scheduler(graph, keys, **kwargs):
        computes.append(keys)
        return dask.get(graph, keys, **kwargs)

    with dask.config.set(scheduler=recording_scheduler):
        found = reader.get_angles()

    times = reader._times_as_np_datetime64[:, np.newaxis]
    tle1, tle2 = reader.get_tle_lines()
    orbit = Orbital(reader.spacecrafts_orbital[reader.spacecraft_id], line1=tle1, line2=tle2)
    sat_azi, sat_elev = orbit.get_observer_look(times, reader.lons, reader.lats, 0)
    sun_azi = np.rad2deg(astronomy.get_alt_az(times, reader.lons, reader.lats)[1])
    expected = (centered_modulus(sat_azi, 360.0), 90 - sat_elev, centered_modulus(sun_azi, 360.0),
                astronomy.sun_zenith_angle(times, reader.lons, reader.lats),
                get_absolute_azimuth_angle_diff(sun_azi, sat_azi))
    assert len(computes) == 1
    for angle, wanted in zip(found, expected, strict=True):
        wanted[reader.mask] = np.nan
        np.testing.assert_allclose(angle, wanted, rtol=0, atol=1e-10)
