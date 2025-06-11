#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 Abhay Devasthale and Martin Raspaud

# Author(s):

#   Abhay Devasthale <abhay.devasthale@smhi.se>
#   Martin Raspaud <martin.raspaud@smhi.se>
#   Adam Dybbroe <adam.dybbroe@smhi.se>

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

"""Reader for GAC KLM data."""

from __future__ import print_function

import logging

import numpy as np

from pygac.gac_reader import GACReader
from pygac.klm_reader import KLMReader, main_klm

LOG = logging.getLogger(__name__)

# video data object

scanline = np.dtype([("scan_line_number", ">u2"),
                     ("scan_line_year", ">u2"),
                     ("scan_line_day_of_year", ">u2"),
                     ("satellite_clock_drift_delta", ">i2"),
                     ("scan_line_utc_time_of_day", ">u4"),
                     ("scan_line_bit_field", ">u2"),
                     ("zero_fill0", ">i2", (5, )),
                     # QUALITY INDICATORS
                     ("quality_indicator_bit_field", ">u4"),
                     ("scan_line_quality_flags", [("reserved", ">u1"),
                                                  ("time_problem_code", ">u1"),
                                                  ("calibration_problem_code", ">u1"),
                                                  ("earth_location_problem_code", ">u1")]),
                     ("calibration_quality_flags", ">u2", (3, )),
                     ("count_of_bit_errors_in_frame_sync", ">u2"),
                     ("zero_fill1", ">i4", (2, )),
                     # CALIBRATION COEFFICIENTS
                     ("visible_operational_cal_ch_1_slope_1", ">i4"),
                     ("visible_operational_cal_ch_1_intercept_1", ">i4"),
                     ("visible_operational_cal_ch_1_slope_2", ">i4"),
                     ("visible_operational_cal_ch_1_intercept_2", ">i4"),
                     ("visible_operational_cal_ch_1_intersection", ">i4"),
                     ("visible_test_cal_ch_1_slope_1", ">i4"),
                     ("visible_test_cal_ch_1_intercept_1", ">i4"),
                     ("visible_test_cal_ch_1_slope_2", ">i4"),
                     ("visible_test_cal_ch_1_intercept_2", ">i4"),
                     ("visible_test_cal_ch_1_intersection", ">i4"),
                     ("visible_prelaunch_cal_ch_1_slope_1", ">i4"),
                     ("visible_prelaunch_cal_ch_1_intercept_1", ">i4"),
                     ("visible_prelaunch_cal_ch_1_slope_2", ">i4"),
                     ("visible_prelaunch_cal_ch_1_intercept_2", ">i4"),
                     ("visible_prelaunch_cal_ch_1_intersection", ">i4"),
                     ("visible_operational_cal_ch_2_slope_1", ">i4"),
                     ("visible_operational_cal_ch_2_intercept_1", ">i4"),
                     ("visible_operational_cal_ch_2_slope_2", ">i4"),
                     ("visible_operational_cal_ch_2_intercept_2", ">i4"),
                     ("visible_operational_cal_ch_2_intersection", ">i4"),
                     ("visible_test_cal_ch_2_slope_1", ">i4"),
                     ("visible_test_cal_ch_2_intercept_1", ">i4"),
                     ("visible_test_cal_ch_2_slope_2", ">i4"),
                     ("visible_test_cal_ch_2_intercept_2", ">i4"),
                     ("visible_test_cal_ch_2_intersection", ">i4"),
                     ("visible_prelaunch_cal_ch_2_slope_1", ">i4"),
                     ("visible_prelaunch_cal_ch_2_intercept_1", ">i4"),
                     ("visible_prelaunch_cal_ch_2_slope_2", ">i4"),
                     ("visible_prelaunch_cal_ch_2_intercept_2", ">i4"),
                     ("visible_prelaunch_cal_ch_2_intersection", ">i4"),
                     ("visible_operational_cal_ch_3a_slope_1", ">i4"),
                     ("visible_operational_cal_ch_3a_intercept_1", ">i4"),
                     ("visible_operational_cal_ch_3a_slope_2", ">i4"),
                     ("visible_operational_cal_ch_3a_intercept_2", ">i4"),
                     ("visible_operational_cal_ch_3a_intersection", ">i4"),
                     ("visible_test_cal_ch_3a_slope_1", ">i4"),
                     ("visible_test_cal_ch_3a_intercept_1", ">i4"),
                     ("visible_test_cal_ch_3a_slope_2", ">i4"),
                     ("visible_test_cal_ch_3a_intercept_2", ">i4"),
                     ("visible_test_cal_ch_3a_intersection", ">i4"),
                     ("visible_prelaunch_cal_ch_3a_slope_1", ">i4"),
                     ("visible_prelaunch_cal_ch_3a_intercept_1", ">i4"),
                     ("visible_prelaunch_cal_ch_3a_slope_2", ">i4"),
                     ("visible_prelaunch_cal_ch_3a_intercept_2", ">i4"),
                     ("visible_prelaunch_cal_ch_3a_intersection", ">i4"),
                     ("ir_operational_cal_ch_3b_coefficient_1", ">i4"),
                     ("ir_operational_cal_ch_3b_coefficient_2", ">i4"),
                     ("ir_operational_cal_ch_3b_coefficient_3", ">i4"),
                     ("ir_test_cal_ch_3b_coefficient_1", ">i4"),
                     ("ir_test_cal_ch_3b_coefficient_2", ">i4"),
                     ("ir_test_cal_ch_3b_coefficient_3", ">i4"),
                     ("ir_operational_cal_ch_4_coefficient_1", ">i4"),
                     ("ir_operational_cal_ch_4_coefficient_2", ">i4"),
                     ("ir_operational_cal_ch_4_coefficient_3", ">i4"),
                     ("ir_test_cal_ch_4_coefficient_1", ">i4"),
                     ("ir_test_cal_ch_4_coefficient_2", ">i4"),
                     ("ir_test_cal_ch_4_coefficient_3", ">i4"),
                     ("ir_operational_cal_ch_5_coefficient_1", ">i4"),
                     ("ir_operational_cal_ch_5_coefficient_2", ">i4"),
                     ("ir_operational_cal_ch_5_coefficient_3", ">i4"),
                     ("ir_test_cal_ch_5_coefficient_1", ">i4"),
                     ("ir_test_cal_ch_5_coefficient_2", ">i4"),
                     ("ir_test_cal_ch_5_coefficient_3", ">i4"),
                     # NAVIGATION
                     ("computed_yaw_steering", ">i2", (3,)),  # only in version 5
                     ("total_applied_attitude_correction",
                      ">i2", (3,)),  # only in version 5
                     ("navigation_status_bit_field", ">u4"),
                     ("time_associated_with_tip_euler_angles", ">u4"),
                     ("tip_euler_angles", ">i2", (3, )),
                     ("spacecraft_altitude_above_reference_ellipsoid", ">u2"),
                     ("angular_relationships", ">i2", (153, )),
                     ("zero_fill3", ">i2", (3, )),
                     ("earth_location", [("lats", ">i4"),
                                         ("lons", ">i4")], (51,)),
                     ("zero_fill4", ">i4", (2, )),
                     # HRPT MINOR FRAME TELEMETRY
                     ("frame_sync", ">u2", (6, )),
                     ("id", ">u2", (2, )),
                     ("time_code", ">u2", (4, )),
                     ("telemetry", [("ramp_calibration", ">u2", (5, )),
                                    ("PRT", ">u2", (3, )),
                                    ("ch3_patch_temp", ">u2"),
                                    ("spare", ">u2"), ]),
                     ("back_scan", ">u2", (30, )),
                     ("space_data", ">u2", (50, )),
                     ("sync_delta", ">u2"),
                     ("zero_fill5", ">i2"),
                     # AVHRR SENSOR DATA
                     ("sensor_data", ">u4", (682, )),
                     ("zero_fill6", ">i4", (2, )),
                     # DIGITAL B TELEMETRY
                     ("invalid_word_bit_flags1", ">u2"),
                     ("avhrr_digital_b_data", ">u2"),
                     ("zero_fill7", ">i4", (3, )),
                     # ANALOG HOUSEKEEPING DATA (TIP)
                     ("invalid_word_bit_flags2", ">u4"),
                     ("patch_temperature_range", ">u1"),
                     ("patch_temperature_extended", ">u1"),
                     ("patch_power", ">u1"),
                     ("radiator_temperature", ">u1"),
                     ("blackbody_temperature_1", ">u1"),
                     ("blackbody_temperature_2", ">u1"),
                     ("blackbody_temperature_3", ">u1"),
                     ("blackbody_temperature_4", ">u1"),
                     ("electronics_current", ">u1"),
                     ("motor_current", ">u1"),
                     ("earth_shield_position", ">u1"),
                     ("electronics_temperature", ">u1"),
                     ("cooler_housing_temperature", ">u1"),
                     ("baseplate_temperature", ">u1"),
                     ("motor_housing_temperature", ">u1"),
                     ("a/d_converter_temperature", ">u1"),
                     ("detector_#4_bias_voltage", ">u1"),
                     ("detector_#5_bias_voltage", ">u1"),
                     ("blackbody_temperature_channel3b", ">u1"),
                     ("blackbody_temperature_channel4", ">u1"),
                     ("blackbody_temperature_channel5", ">u1"),
                     ("reference_voltage", ">u1"),
                     ("zero_fill8", ">i2", (3, )),
                     # CLOUDS FROM AVHRR (CLAVR)
                     ("reserved0", ">u4"),
                     ("reserved1", ">u4"),
                     ("reserved2", ">u2", (52, )),
                     # FILLER
                     ("zero_fill9", ">i4", (112, ))])


class GACKLMReader(GACReader, KLMReader):
    """The GAC KLM reader class.

    The offset attribute tells where in the file the scanline data starts.
    """
    lonlat_sample_points = np.arange(4.5, 405, 8)

    def __init__(self, *args, **kwargs):
        """Init the GAC KLM reader."""
        GACReader.__init__(self, *args, **kwargs)
        self.scanline_type = scanline
        self.offset = 4608


def main(filename, start_line, end_line):
    """Generate a l1c file."""
    return main_klm(GACKLMReader, filename, start_line, end_line)


if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2], sys.argv[3])
