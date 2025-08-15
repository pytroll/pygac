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
from pygac.reader import NoTLEData


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
        """Get the telemetry."""
        prt = 51 * np.ones(self.along_track)  # prt threshold is 50
        prt[::5] = 0
        ict = 101 * np.ones((self.along_track, 3))  # ict threshold is 100
        space = 101 * np.ones((self.along_track, 3))  # space threshold is 100
        total_space = 101 * np.ones((self.along_track, 10, 3))  # space threshold is 100
        total_ict = 101 * np.ones((self.along_track, 10, 3))  # space threshold is 100
        return prt, ict, space, total_space, total_ict

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
        """Get the telemetry."""
        prt = 51 * np.ones(self.along_track)  # prt threshold is 50
        prt[::5] = 0
        ict = 101 * np.ones((self.along_track, 3))  # ict threshold is 100
        space = 101 * np.ones((self.along_track, 3))  # space threshold is 100
        total_space = 101 * np.ones((self.along_track, 10, 3))  # space threshold is 100
        total_ict = 101 * np.ones((self.along_track, 10, 3))  # space threshold is 100
        return prt, ict, space, total_space, total_ict

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
            self.reader._times_as_np_datetime64 = utcs
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
        self.reader._times_as_np_datetime64 = np.array(
            [315748035469, 315748359969,
             315751135469, 315754371969,
             315754371969]).astype("datetime64[ms]")
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
        self.reader._times_as_np_datetime64 = np.array([315748035469, 315748359969,
                                     315751135469, 315754371969,
                                     315754371969]).astype("datetime64[ms]")
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
    assert dataset["prt_counts"].shape == (3,)
    assert dataset["ict_counts"].shape == (3, 3)
    assert dataset["space_counts"].shape == (3, 3)
#
# Code seems to require interpolation to with with the new addition
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
    assert lons[0, 0] == pytest.approx(-4.756486)
    assert lons[0, -1] == pytest.approx(79.505688)


def test_recomputing_lonlats_with_time_offset(pod_file_with_tbm_header, pod_tle):
    """Test computing lons and lats from TLE data with time offset."""
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name, compute_lonlats_from_tles=True)
    reader.read(pod_file_with_tbm_header)
    lons, lats = reader._compute_lonlats(time_offset=np.timedelta64(500, "ms"))
    assert lons[0, 0] == pytest.approx(-4.714710005688344)
    assert lons[0, -1] == pytest.approx(79.44587709203307)


def test_georeferencing_with_first_guess(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """Test getting a first guess for georeferencing from the file lon/lats."""
    expected_time_offset = 2  # seconds
    def mock_cal(counts, *args, **kwargs):
        return counts * 1.0
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal", mock_cal)

    def mock_disp(*args):
        return 0, (0, 0, 0), ([10000], [1000])
    from georeferencer import georeferencer
    monkeypatch.setattr(georeferencer, "get_swath_displacement", mock_disp)
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
    assert tdiff < np.timedelta64(2, "ms")


def test_georeferencing_fails(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """Test georeferencing."""

    def mock_cal(counts, *args, **kwargs):
        return counts * 1.0
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal", mock_cal)

    def mock_disp(*args):
        return 0, (0, 0, 0), ([10000], [10000])
    from georeferencer import georeferencer
    monkeypatch.setattr(georeferencer, "get_swath_displacement", mock_disp)
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name, compute_lonlats_from_tles=True,
                          reference_image="some_world_image.tif")
    reader.read(pod_file_with_tbm_header)
    with pytest.raises(RuntimeError):
        _ = reader.get_calibrated_dataset()


def test_georeferencing(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """Test georeferencing."""

    def mock_cal(counts, *args, **kwargs):
        return counts * 1.0
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal", mock_cal)

    def mock_disp(*args):
        return 0.5, (0, 0, 0), ([10000], [1000])
    from georeferencer import georeferencer
    monkeypatch.setattr(georeferencer, "get_swath_displacement", mock_disp)
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name, compute_lonlats_from_tles=True,
                          reference_image="some_world_image.tif")
    reader.read(pod_file_with_tbm_header)
    dataset = reader.get_calibrated_dataset()
    assert dataset.attrs["max_scan_angle"] == 55.37
    lons = dataset["longitude"].values
    assert lons[0, 0] == pytest.approx(-4.714710005688344)
    assert lons[0, -1] == pytest.approx(79.44587709203307)

def test_orthocorrection(pod_file_with_tbm_header, pod_tle, monkeypatch):
    """Test computing lons and lats from TLE data."""

    def mock_cal(counts, *args, **kwargs):
        return counts * 1.0
    import pygac.calibration.noaa
    monkeypatch.setattr(pygac.calibration.noaa, "calibrate_thermal", mock_cal)

    def mock_disp(*args):
        return 0.5, (0, 0, 0), ([10000], [1000])
    from georeferencer import georeferencer
    monkeypatch.setattr(georeferencer, "get_swath_displacement", mock_disp)
    reader = LACPODReader(tle_dir=pod_tle.parent, tle_name=pod_tle.name, compute_lonlats_from_tles=True,
                          reference_image="some_world_image.tif", dem="dem_file.tif")
    reader.read(pod_file_with_tbm_header)


    def mock_orthocorrection(calibrated_ds, *args):
        calibrated_ds["tc_lons"] = calibrated_ds["longitude"] + 0.0001
        calibrated_ds["tc_lats"] = calibrated_ds["latitude"] + 0.0001
        return calibrated_ds

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
