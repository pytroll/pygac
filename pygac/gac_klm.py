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

"""Read a gac file.
Reads L1b GAC data from KLM series of satellites (NOAA-15 and later) and does most of the computations.
Format specification can be found here:
http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/c8/sec83142-1.htm

"""

import numpy as np
import geotiepoints as gtp
import datetime
from pyorbital import astronomy
import pygac.calibrate_klm as cal_klm
import gac_io
from pyorbital.orbital import Orbital

import logging
LOG = logging.getLogger(__name__)

import ConfigParser
import os

try:
    CONFIG_FILE = os.environ['PYGAC_CONFIG_FILE']
except KeyError:
    LOG.exception('Environment variable PYGAC_CONFIG_FILE not set!')
    raise

if not os.path.exists(CONFIG_FILE) or not os.path.isfile(CONFIG_FILE):
    raise IOError(str(CONFIG_FILE) + " pointed to by the environment " + 
                  "variable PYGAC_CONFIG_FILE is not a file or does not exist!")


# setting up constants used in the processing

AVHRR_SWATH_WIDTH_DEG=55.385
sat_altitude=850000.0
M_PI=3.14159265359
AHA_DEG_TO_RAD=M_PI/180.
AHA_RAD_TO_DEG=180./M_PI
earth_radius=6370997.0
wholeday_ms=24.0*3600.0*1000.0
GAC_NEAR_NADIR_POSITION=207
MISSING_DATA = -32001


# GAC header object

header = np.dtype([("data_set_creation_site_id", "S3"),
                   ("ascii_blank_=_x20", "S1"),
                   ("noaa_level_1b_format_version_number:", ">u2"),
                   ("noaa_level_1b_format_version_year", ">u2"),
                   ("noaa_level_1b_format_version_day_of_year", ">u2"),
                   ("reserved_for_logical_record_length", ">u2"),
                   ("reserved_for_block_size", ">u2"),
                   ("count_of_header_records", ">u2"),
                   ("zero_fill0", ">i2", (3, )),
                   ("data_set_name", "S42"),
                   ("processing_block_identification", "S8"),
                   ("noaa_spacecraft_identification_code", ">u2"),
                   ("instrument_id", ">u2"),
                   ("data_type_code", ">u2"),
                   ("tip_source_code", ">u2"),
                   ("start_of_data_set_day_count_starting_from_0_at_00h,_1_jan_1950", ">u4"),
                   ("start_of_data_set_year", ">u2"),
                   ("start_of_data_set_day_of_year", ">u2"),
                   ("start_of_data_set_utc_time_of_day", ">u4"),
                   ("end_of_data_set_day_count_starting_from_0_at_00h,_1_jan_1950", ">u4"),
                   ("end_of_data_set_year", ">u2"),
                   ("end_of_data_set_day_of_year", ">u2"),
                   ("end_of_data_set_utc_time_of_day", ">u4"),
                   ("year_of_last_cpids_update", ">u2"),
                   ("day_of_year_of_last_cpids_update", ">u2"),
                   ("zero_fill1", ">i2", (4, )),
# DATA SET QUALITY INDICATORS
                   ("instrument_status", ">u4"),
                   ("zero_fill2", ">i2"),
                   ("record_number_of_status_change", ">u2"),
                   ("second_instrument_status", ">u4"),
                   ("count_of_data_records", ">u2"),
                   ("count_of_calibrated,_earth_located_scan_lines", ">u2"),
                   ("count_of_missing_scan_lines", ">u2"),
                   ("count_of_data_gaps", ">u2"),
                   ("filler0", ">u2"),
                   ("count_of_pacs_detected_tip_parity_errors", ">u2"),
                   ("filler1", ">u2"),
                   ("time_sequence_error", ">u2"),
                   ("time_sequence_error_code", ">u2"),
                   ("socc_clock_update_indicator", ">u2"),
                   ("earth_location_error_indicator", ">u2"),
                   ("earth_location_error_code", ">u2"),
                   ("pacs_status_bit_field", ">u2"),
                   ("data_source", ">u2"),
                   ("zero_fill3", ">i4"),
                   ("reserved_for_the_ingester", "S8"),
                   ("reserved_for_decommutation", "S8"),
                   ("zero_fill4", ">i2", (5, )),
# CALIBRATION
                   ("ramp_auto_calibration_indicators_bit_field", ">u2"),
                   ("year_of_most_recent_solar_channel_calibration", ">u2"),
                   ("day_of_year_of_most_recent_solar_channel_calibration", ">u2"),
                   ("primary_calibration_algorithm_id", ">u2"),
                   ("primary_calibration_algorithm_selected_options", ">u2"),
                   ("secondary_calibration_algorithm_id", ">u2"),
                   ("secondary_calibration_algorithm_selected_options", ">u2"),
                   ("ir_target_temperature_1_conversion_coefficient_1", ">i2"),
                   ("ir_target_temperature_1_conversion_coefficient_2", ">i2"),
                   ("ir_target_temperature_1_conversion_coefficient_3", ">i2"),
                   ("ir_target_temperature_1_conversion_coefficient_4", ">i2"),
                   ("ir_target_temperature_1_conversion_coefficient_5", ">i2"),
                   ("ir_target_temperature_1_conversion_coefficient_6", ">i2"),
                   ("ir_target_temperature_2_conversion_coefficient_1", ">i2"),
                   ("ir_target_temperature_2_conversion_coefficient_2", ">i2"),
                   ("ir_target_temperature_2_conversion_coefficient_3", ">i2"),
                   ("ir_target_temperature_2_conversion_coefficient_4", ">i2"),
                   ("ir_target_temperature_2_conversion_coefficient_5", ">i2"),
                   ("ir_target_temperature_2_conversion_coefficient_6", ">i2"),
                   ("ir_target_temperature_3_conversion_coefficient_1", ">i2"),
                   ("ir_target_temperature_3_conversion_coefficient_2", ">i2"),
                   ("ir_target_temperature_3_conversion_coefficient_3", ">i2"),
                   ("ir_target_temperature_3_conversion_coefficient_4", ">i2"),
                   ("ir_target_temperature_3_conversion_coefficient_5", ">i2"),
                   ("ir_target_temperature_3_conversion_coefficient_6", ">i2"),
                   ("ir_target_temperature_4_conversion_coefficient_1", ">i2"),
                   ("ir_target_temperature_4_conversion_coefficient_2", ">i2"),
                   ("ir_target_temperature_4_conversion_coefficient_3", ">i2"),
                   ("ir_target_temperature_4_conversion_coefficient_4", ">i2"),
                   ("ir_target_temperature_4_conversion_coefficient_5", ">i2"),
                   ("ir_target_temperature_4_conversion_coefficient_6", ">i2"),
                   ("zero_fill5", ">i4", (2, )),
# RADIANCE CONVERSION
                   ("ch_1_solar_filtered_irradiance_in_wavelength", ">i4"),
                   ("ch_1_equivalent_filter_width_in_wavelength", ">i4"),
                   ("ch_2_solar_filtered_irradiance_in_wavelength", ">i4"),
                   ("ch_2_equivalent_filter_width_in_wavelength", ">i4"),
                   ("ch_3a_solar_filtered_irradiance_in_wavelength", ">i4"),
                   ("ch_3a_equivalent_filter_width_in_wavelength", ">i4"),
                   ("ch_3b_central_wavenumber_(see", ">i4"),
                   ("ch_3b_constant_1", ">i4"),
                   ("ch_3b_constant_2", ">i4"),
                   ("ch_4_central_wavenumber", ">i4"),
                   ("ch_4_constant_1", ">i4"),
                   ("ch_4_constant_2", ">i4"),
                   ("ch_5_central_wavenumber", ">i4"),
                   ("ch_5_constant_1", ">i4"),
                   ("ch_5_constant_2", ">i4"),
                   ("zero_fill6", ">i4", (3, )),
# NAVIGATION
                   ("reference_ellipsoid_model_id", "S8"),
                   ("nadir_earth_location_tolerance", ">u2"),
                   ("earth_location_bit_field", ">u2"),
                   ("zero_fill7", ">i2"),
                   ("constant_roll_attitude_error", ">i2"),
                   ("constant_pitch_attitude_error", ">i2"),
                   ("constant_yaw_attitude_error", ">i2"),
                   ("epoch_year_for_orbit_vector", ">u2"),
                   ("day_of_epoch_year_for_orbit_vector", ">u2"),
                   ("epoch_utc_time_of_day_for_orbit_vector", ">u4"),
                   ("semi-major_axis", ">i4"),
                   ("eccentricity", ">i4"),
                   ("inclination", ">i4"),
                   ("argument_of_perigee", ">i4"),
                   ("right_ascension_of_the_ascending_node", ">i4"),
                   ("mean_anomaly", ">i4"),
                   ("position_vector_x_component", ">i4"),
                   ("position_vector_y_component", ">i4"),
                   ("position_vector_z_component", ">i4"),
                   ("velocity_vector_x-dot_component", ">i4"),
                   ("velocity_vector_y-dot_component", ">i4"),
                   ("velocity_vector_z-dot_component", ">i4"),
                   ("earth_sun_distance_ratio", ">u4"),
                   ("zero_fill8", ">i4", (4, )),
# ANALOG TELEMETRY CONVERSION
                   ("patch_temperature_conversion_coefficient_1", ">i2"),
                   ("patch_temperature_conversion_coefficient_2", ">i2"),
                   ("patch_temperature_conversion_coefficient_3", ">i2"),
                   ("patch_temperature_conversion_coefficient_4", ">i2"),
                   ("patch_temperature_conversion_coefficient_5", ">i2"),
                   ("reserved0", ">i2"),
                   ("patch_temperature_extended_conversion_coefficient_1", ">i2"),
                   ("patch_temperature_extended_conversion_coefficient_2", ">i2"),
                   ("patch_temperature_extended_conversion_coefficient_3", ">i2"),
                   ("patch_temperature_extended_conversion_coefficient_4", ">i2"),
                   ("patch_temperature_extended_conversion_coefficient_5", ">i2"),
                   ("reserved1", ">i2"),
                   ("patch_power_conversion_coefficient_1", ">i2"),
                   ("patch_power_conversion_coefficient_2", ">i2"),
                   ("patch_power_conversion_coefficient_3", ">i2"),
                   ("patch_power_conversion_coefficient_4", ">i2"),
                   ("patch_power_conversion_coefficient_5", ">i2"),
                   ("reserved2", ">i2"),
                   ("radiator_temperature_conversion_coefficient_1", ">i2"),
                   ("radiator_temperature_conversion_coefficient_2", ">i2"),
                   ("radiator_temperature_conversion_coefficient_3", ">i2"),
                   ("radiator_temperature_conversion_coefficient_4", ">i2"),
                   ("radiator_temperature_conversion_coefficient_5", ">i2"),
                   ("reserved3", ">i2"),
                   ("blackbody_temperature_1_conversion_coefficient_1", ">i2"),
                   ("blackbody_temperature_1_conversion_coefficient_2", ">i2"),
                   ("blackbody_temperature_1_conversion_coefficient_3", ">i2"),
                   ("blackbody_temperature_1_conversion_coefficient_4", ">i2"),
                   ("blackbody_temperature_1_conversion_coefficient_5", ">i2"),
                   ("reserved4", ">i2"),
                   ("blackbody_temperature_2_conversion_coefficient_1", ">i2"),
                   ("blackbody_temperature_2_conversion_coefficient_2", ">i2"),
                   ("blackbody_temperature_2_conversion_coefficient_3", ">i2"),
                   ("blackbody_temperature_2_conversion_coefficient_4", ">i2"),
                   ("blackbody_temperature_2_conversion_coefficient_5", ">i2"),
                   ("reserved5", ">i2"),
                   ("blackbody_temperature_3_conversion_coefficient_1", ">i2"),
                   ("blackbody_temperature_3_conversion_coefficient_2", ">i2"),
                   ("blackbody_temperature_3_conversion_coefficient_3", ">i2"),
                   ("blackbody_temperature_3_conversion_coefficient_4", ">i2"),
                   ("blackbody_temperature_3_conversion_coefficient_5", ">i2"),
                   ("reserved6", ">i2"),
                   ("blackbody_temperature_4_conversion_coefficient_1", ">i2"),
                   ("blackbody_temperature_4_conversion_coefficient_2", ">i2"),
                   ("blackbody_temperature_4_conversion_coefficient_3", ">i2"),
                   ("blackbody_temperature_4_conversion_coefficient_4", ">i2"),
                   ("blackbody_temperature_4_conversion_coefficient_5", ">i2"),
                   ("reserved7", ">i2"),
                   ("electronics_current_conversion_coefficient_1", ">i2"),
                   ("electronics_current_conversion_coefficient_2", ">i2"),
                   ("electronics_current_conversion_coefficient_3", ">i2"),
                   ("electronics_current_conversion_coefficient_4", ">i2"),
                   ("electronics_current_conversion_coefficient_5", ">i2"),
                   ("reserved8", ">i2"),
                   ("motor_current_conversion_coefficient_1", ">i2"),
                   ("motor_current_conversion_coefficient_2", ">i2"),
                   ("motor_current_conversion_coefficient_3", ">i2"),
                   ("motor_current_conversion_coefficient_4", ">i2"),
                   ("motor_current_conversion_coefficient_5", ">i2"),
                   ("reserved9", ">i2"),
                   ("earth_shield_position_conversion_coefficient_1", ">i2"),
                   ("earth_shield_position_conversion_coefficient_2", ">i2"),
                   ("earth_shield_position_conversion_coefficient_3", ">i2"),
                   ("earth_shield_position_conversion_coefficient_4", ">i2"),
                   ("earth_shield_position_conversion_coefficient_5", ">i2"),
                   ("reserved10", ">i2"),
                   ("electronics_temperature_conversion_coefficient_1", ">i2"),
                   ("electronics_temperature_conversion_coefficient_2", ">i2"),
                   ("electronics_temperature_conversion_coefficient_3", ">i2"),
                   ("electronics_temperature_conversion_coefficient_4", ">i2"),
                   ("electronics_temperature_conversion_coefficient_5", ">i2"),
                   ("reserved11", ">i2"),
                   ("cooler_housing_temperature_conversion_coefficient_1", ">i2"),
                   ("cooler_housing_temperature_conversion_coefficient_2", ">i2"),
                   ("cooler_housing_temperature_conversion_coefficient_3", ">i2"),
                   ("cooler_housing_temperature_conversion_coefficient_4", ">i2"),
                   ("cooler_housing_temperature_conversion_coefficient_5", ">i2"),
                   ("reserved12", ">i2"),
                   ("baseplate_temperature_conversion_coefficient_1", ">i2"),
                   ("baseplate_temperature_conversion_coefficient_2", ">i2"),
                   ("baseplate_temperature_conversion_coefficient_3", ">i2"),
                   ("baseplate_temperature_conversion_coefficient_4", ">i2"),
                   ("baseplate_temperature_conversion_coefficient_5", ">i2"),
                   ("reserved13", ">i2"),
                   ("motor_housing_temperature_conversion_coefficient_1", ">i2"),
                   ("motor_housing_temperature_conversion_coefficient_2", ">i2"),
                   ("motor_housing_temperature_conversion_coefficient_3", ">i2"),
                   ("motor_housing_temperature_conversion_coefficient_4", ">i2"),
                   ("motor_housing_temperature_conversion_coefficient_5", ">i2"),
                   ("reserved14", ">i2"),
                   ("a/d_converter_temperature_conversion_coefficient_1", ">i2"),
                   ("a/d_converter_temperature_conversion_coefficient_2", ">i2"),
                   ("a/d_converter_temperature_conversion_coefficient_3", ">i2"),
                   ("a/d_converter_temperature_conversion_coefficient_4", ">i2"),
                   ("a/d_converter_temperature_conversion_coefficient_5", ">i2"),
                   ("reserved15", ">i2"),
                   ("detector_#4_bias_voltage_conversion_coefficient_1", ">i2"),
                   ("detector_#4_bias_voltage_conversion_coefficient_2", ">i2"),
                   ("detector_#4_bias_voltage_conversion_coefficient_3", ">i2"),
                   ("detector_#4_bias_voltage_conversion_coefficient_4", ">i2"),
                   ("detector_#4_bias_voltage_conversion_coefficient_5", ">i2"),
                   ("reserved16", ">i2"),
                   ("detector_#5_bias_voltage_conversion_coefficient_1", ">i2"),
                   ("detector_#5_bias_voltage_conversion_coefficient_2", ">i2"),
                   ("detector_#5_bias_voltage_conversion_coefficient_3", ">i2"),
                   ("detector_#5_bias_voltage_conversion_coefficient_4", ">i2"),
                   ("detector_#5_bias_voltage_conversion_coefficient_5", ">i2"),
                   ("reserved17", ">i2"),
                   ("channel_3b_blackbody_view_conversion_coefficient_1", ">i2"),
                   ("channel_3b_blackbody_view_conversion_coefficient_2", ">i2"),
                   ("channel_3b_blackbody_view_conversion_coefficient_3", ">i2"),
                   ("channel_3b_blackbody_view_conversion_coefficient_4", ">i2"),
                   ("channel_3b_blackbody_view_conversion_coefficient_5", ">i2"),
                   ("reserved18", ">i2"),
                   ("channel_4_blackbody_view_conversion_coefficient_1", ">i2"),
                   ("channel_4_blackbody_view_conversion_coefficient_2", ">i2"),
                   ("channel_4_blackbody_view_conversion_coefficient_3", ">i2"),
                   ("channel_4_blackbody_view_conversion_coefficient_4", ">i2"),
                   ("channel_4_blackbody_view_conversion_coefficient_5", ">i2"),
                   ("reserved19", ">i2"),
                   ("channel_5_blackbody_view_conversion_coefficient_1", ">i2"),
                   ("channel_5_blackbody_view_conversion_coefficient_2", ">i2"),
                   ("channel_5_blackbody_view_conversion_coefficient_3", ">i2"),
                   ("channel_5_blackbody_view_conversion_coefficient_4", ">i2"),
                   ("channel_5_blackbody_view_conversion_coefficient_5", ">i2"),
                   ("reserved20", ">i2"),
                   ("reference_voltage_conversion_coefficient_1", ">i2"),
                   ("reference_voltage_conversion_coefficient_2", ">i2"),
                   ("reference_voltage_conversion_coefficient_3", ">i2"),
                   ("reference_voltage_conversion_coefficient_4", ">i2"),
                   ("reference_voltage_conversion_coefficient_5", ">i2"),
                   ("reserved21", ">i2")])
# FILLER
#                  ("zero_fill9", ">i4", (4608 - 688, ))])


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
                     ("scan_line_quality_flags", ">u4"),
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
                     ("zero_fill2", ">i4", (3, )),
                     # NAVIGATION
                     ("navigation_status_bit_field", ">u4"),
                     ("time_associated_with_tip_euler_angles", ">u4"),
                     ("tip_euler_angles", ">i2", (3, )),
                     ("spacecraft_altitude_above_reference_ellipsoid", ">u2"),
                     ("angular_relationships", ">i2", (153, )),
                     ("zero_fill3", ">i2", (3, )),
                     ("earth_location", ">i4", (102, )),
                     ("zero_fill4", ">i4", (2, )),
                     # HRPT MINOR FRAME TELEMETRY
                     ("frame_sync", ">u2", (6, )),
                     ("id", ">u2", (2, )),
                     ("time_code", ">u2", (4, )),
                     ("telemetry", ">u2", (10, )),
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
                     ("word_1_patch_temperature", ">u1", (22, )),
                     ("zero_fill8", ">i2", (3, )),
                     # CLOUDS FROM AVHRR (CLAVR)
                     ("reserved0", ">u4"),
                     ("reserved1", ">u4"),
                     ("reserved2", ">u2", (52, )),
                     # FILLER
                     ("zero_fill9", ">i4", (112, ))])



def main(filename):
    # reading raw L1b data

    with open(filename) as fd_:
        head = np.fromfile(fd_, dtype=header, count=1)
        fd_.seek(4608, 0)
        scans = np.fromfile(fd_, dtype=scanline)

    spacecraft_id=int(head["noaa_spacecraft_identification_code"])
    spacecrafts = {4: 'noaa15',
                   2: 'noaa16',
		   6: 'noaa17',
		   7: 'noaa18',
		   8: 'noaa19',
		   12: 'metop02',
                   }
    
    try:
        satellite_name = spacecrafts[spacecraft_id]
    except KeyError:
        raise KeyError("Wrong satellite id: " + str(spacecraft_id))


    conf = ConfigParser.ConfigParser()
    try:
        conf.read(CONFIG_FILE)
    except ConfigParser.NoSectionError:
        LOG.exception('Failed reading configuration file: ' + str(CONFIG_FILE))
        raise

    values = {"satname": satellite_name, }

    options = {}
    for option, value in conf.items('tle', raw = True):
        options[option] = value

    tle_filename = os.path.join(options['tledir'], 
                                options["tlename"] % values)
    LOG.info('TLE filename = ' + str(tle_filename))


    # unpacking raw data to get counts

    packed_data = scans["sensor_data"]
    gac_counts = np.zeros((int(head["count_of_data_records"]),
                     409 * 5))
    gac_counts[:, 0::3] = (packed_data & (1023 << 20)) >> 20
    gac_counts[:, 1::3] = (packed_data & (1023 << 10)) >> 10
    gac_counts[:, 2::3] = (packed_data & 1023)[:, :-1]



    # interpolating lat-on points using PYTROLL geotiepoints

    arrLat = np.zeros((float(head["count_of_data_records"]), 51))
    arrLon = np.zeros((float(head["count_of_data_records"]), 51))
    arrLat = scans["earth_location"][:, 0::2] / pow(10.0, 4.0)
    arrLon = scans["earth_location"][:, 1::2] / pow(10.0, 4.0)

    arrLon_full, arrLat_full = gtp.Gac_Lat_Lon_Interpolator(arrLon,arrLat)


    # getting time information and calculating solar zenith angle for the entire orbit

    arrYear = scans["scan_line_year"]
    arrJDay = scans["scan_line_day_of_year"]
    arrUTC = scans["scan_line_utc_time_of_day"]

    arrSZA = np.zeros((int(head["count_of_data_records"]), 409))
    for i in range(head["count_of_data_records"]):
            temp_utc = (datetime.datetime(arrYear[i], 1, 1) 
                        + datetime.timedelta(int(arrJDay[i]) - 1) 
                        + datetime.timedelta(milliseconds=int(arrUTC[i])))
            arrSZA[i,:]= astronomy.sun_zenith_angle(temp_utc, arrLon_full[i,:], arrLat_full[i,:])




    # calculating satellite zenith angle

    pixel_pos = np.zeros((int(head["count_of_data_records"]), 409))
    for i in range(head["count_of_data_records"]):
            pixel_pos[i,:] = np.arange(0,409,1)
    #pixel_pos = np.tile(np.arange(409), (head["count_of_data_records"], 1))

    scan_angle = AVHRR_SWATH_WIDTH_DEG * np.divide(np.absolute(pixel_pos - 205), 205.0)
    #scan_angle = AVHRR_SWATH_WIDTH_DEG * abs(pixel_pos - 205) / 205.0
    # this isn't needed
    arrSTZ = np.zeros((int(head["count_of_data_records"]), 409))
    arrSTZ = np.arcsin((1.0+sat_altitude/earth_radius)*np.sin(scan_angle*AHA_DEG_TO_RAD))*AHA_RAD_TO_DEG;



    # calculating solar azimuth angle

    arrSAA = np.zeros((int(head["count_of_data_records"]), 409))
    for i in range(head["count_of_data_records"]):
            temp_utc = datetime.datetime(arrYear[i], 1, 1) + datetime.timedelta(int(arrJDay[i]) - 1) + datetime.timedelta(milliseconds=int(arrUTC[i]));
            temp_alt, temp_suna = astronomy.get_alt_az(temp_utc, arrLon_full[i,:], arrLat_full[i,:]);	
            temp_suna=temp_suna*180.0/M_PI
            suna=temp_suna*0.0
            ii=np.where(temp_suna<0.0)
            suna[ii]=temp_suna[ii]+180.0
            jj=np.where(temp_suna>=0.0)
            suna[jj]=temp_suna[jj]-180.0
            arrSAA[i,:]=suna


    # calculating satellite azimuth angle

    tle_name = tle_filename # './noaa_tle/TLE_%s.txt' % satellite_name
    if arrYear[1]<=1999:
	yearJday = '%02d%03d.' % (arrYear[1]%1900,arrJDay[1])
    else:
	yearJday = '%02d%03d.' % (arrYear[1]%2000,arrJDay[1])

    tlep=open(tle_name,'rt')
    tle_data=tlep.readlines()
    tlep.close()

    def find_tle_index(tle_data, yearJday):
    	for iindex, stringYearJday in enumerate(tle_data):
        	if yearJday in stringYearJday:
              		return iindex
    	return -1
   
    try:
    	iindex=find_tle_index(tle_data,yearJday)
    	tle1=tle_data[iindex]
    	tle2=tle_data[iindex+1]
    except KeyError:
    	print "Could't find TLE - exit";
        sys.exit(1);
    
   
    
    spacecraft_id=int(head["noaa_spacecraft_identification_code"])
    spacecrafts_orbital = {4: 'noaa 15',
                   2: 'noaa 16',
		   6: 'noaa 17',
		   7: 'noaa 18',
		   8: 'noaa 19',
		   12: 'metop 02',
                   }
    
    try:
        satellite_name_orbital = spacecrafts_orbital[spacecraft_id]
    except KeyError:
        print "wrong satellite id - exit";
        sys.exit(1);
    
    orb = Orbital(satellite_name_orbital, line1=tle1, line2=tle2)  

    arrSTA = np.zeros((int(head["count_of_data_records"]), 409), dtype='float')
    for i in range(head["count_of_data_records"]):
            temp_utc = datetime.datetime(int(arrYear[i]), 1, 1) + datetime.timedelta(int(arrJDay[i]) - 1) + datetime.timedelta(milliseconds=int(arrUTC[i]));

    	    arr_psta, arr_ele = orb.get_observer_look(temp_utc, arrLon_full[i,:], arrLat_full[i,:], 0)
    	    arrSTA[i,:]=arr_psta;

    
    arrSTA[:,205]=arrSTA[:,204]


    # calculating relative azimuth angle

    arrRAA = np.zeros((int(head["count_of_data_records"]), 409))
    arrRAA = np.absolute(arrSTA-arrSAA)	
    arrRAA = np.where(arrRAA>180.0, 360.0-arrRAA, arrRAA)


    # scaling angles and lat-lon values

    arrSTA=arrSTA-180.0;
    arrRAA=np.where(arrRAA<0.0,-1.0*arrRAA,arrRAA);
    arrRAA=180.0-arrRAA;

    arrSZA=arrSZA*100.0
    arrSTZ=arrSTZ*100.0
    arrSAA=arrSAA*100.0
    arrSTA=arrSTA*100.0
    arrRAA=arrRAA*100.0

    arrLat_full=arrLat_full*100.0
    arrLon_full=arrLon_full*100.0


    # Earth-Sun distance correction factor 
    corr = 1.0 - 0.0334*np.cos(2.0*M_PI*(head["start_of_data_set_day_of_year"]-2)/365.25)


    # Calibrating solar channels

    channel3_switch=np.zeros(int(head["count_of_data_records"]))
    channel3_switch=scans["scan_line_bit_field"][:] & 3
    ref1,ref2,ref3=cal_klm.calibrate_solar(gac_counts, int(head["start_of_data_set_year"]), int(head["start_of_data_set_day_of_year"]), int(head["noaa_spacecraft_identification_code"]), channel3_switch, corr, int(head["count_of_data_records"]))

    ii=np.where(channel3_switch==0)
    ref3[ii]=MISSING_DATA

    # Calibrating thermal channels

    # getting PRT counts

    prt_counts=np.zeros(int(head["count_of_data_records"]))
    prt_counts=np.mean(scans["telemetry"][:,5:8], axis=1)
    

    # getting ICT counts

    ict_counts=np.zeros((int(head["count_of_data_records"]),3))
    ict_counts[:,0]=np.mean(scans["back_scan"][:,0::3],axis=1)
    ict_counts[:,1]=np.mean(scans["back_scan"][:,1::3],axis=1)
    ict_counts[:,2]=np.mean(scans["back_scan"][:,2::3],axis=1)

    # getting space counts

    space_counts=np.zeros((int(head["count_of_data_records"]),3))
    space_counts[:,0]=np.mean(scans["space_data"][:,2::5],axis=1)
    space_counts[:,1]=np.mean(scans["space_data"][:,3::5],axis=1)
    space_counts[:,2]=np.mean(scans["space_data"][:,4::5],axis=1)


    # calibrating channels 3b, 4 and 5

    bt3=cal_klm.calibrate_thermal(gac_counts[:,2::5], prt_counts, ict_counts[:,0], space_counts[:,0], int(head["count_of_data_records"]), int(head["noaa_spacecraft_identification_code"]), channel=3)
    bt4=cal_klm.calibrate_thermal(gac_counts[:,3::5], prt_counts, ict_counts[:,1], space_counts[:,1], int(head["count_of_data_records"]), int(head["noaa_spacecraft_identification_code"]), channel=4)
    bt5=cal_klm.calibrate_thermal(gac_counts[:,4::5], prt_counts, ict_counts[:,2], space_counts[:,2], int(head["count_of_data_records"]), int(head["noaa_spacecraft_identification_code"]), channel=5)

    bt3=(bt3-273.15)*100.0
    bt4=(bt4-273.15)*100.0
    bt5=(bt5-273.15)*100.0
    ii=np.where(channel3_switch==1)
    bt3[ii]=MISSING_DATA
    
    # masking out below-quality scanlines

    scanline_quality=np.zeros((int(head["count_of_data_records"]),3))
    scanline_quality[:,0]=scans["quality_indicator_bit_field"]>>31
    scanline_quality[:,1]=(scans["quality_indicator_bit_field"]<<3)>>31
    scanline_quality[:,2]=(scans["quality_indicator_bit_field"]<<4)>>31

    ii=np.where((scanline_quality[:,0]==1) | (scanline_quality[:,1]==1) | (scanline_quality[:,2]==1))
    arrLat_full[ii]=MISSING_DATA
    arrLon_full[ii]=MISSING_DATA
    ref1[ii]=MISSING_DATA
    ref2[ii]=MISSING_DATA
    ref3[ii]=MISSING_DATA
    bt3[ii]=MISSING_DATA
    bt4[ii]=MISSING_DATA
    bt5[ii]=MISSING_DATA
    arrSZA[ii]=MISSING_DATA
    arrSTZ[ii]=MISSING_DATA
    arrSAA[ii]=MISSING_DATA
    arrSTA[ii]=MISSING_DATA
    arrRAA[ii]=MISSING_DATA
 
    ii=np.where(channel3_switch==2)
    bt3[ii]=MISSING_DATA
    ref3[ii]=MISSING_DATA

    ii=np.where(ref1<0);   
    ref1[ii]=MISSING_DATA;
    ii=np.where(ref2<0); 
    ref2[ii]=MISSING_DATA;

    # writing out calibrated AVHRR channel data and various sun-sat angles        


    t=datetime.datetime(int(head["start_of_data_set_year"]), 1, 1) + datetime.timedelta(int(head["start_of_data_set_day_of_year"]) - 1) + datetime.timedelta(milliseconds=int(head["start_of_data_set_utc_time_of_day"]))
    tenth_s = int(t.microsecond/100000)
    startdate = '%d%02d%02d' % (t.year,t.month,t.day)
    starttime = '%02d%02d%02d%01d' % (t.hour,t.minute,t.second,tenth_s)

    t = datetime.datetime(int(head["end_of_data_set_year"]), 1, 1) + datetime.timedelta(int(head["end_of_data_set_day_of_year"]) - 1) + datetime.timedelta(milliseconds=int(head["end_of_data_set_utc_time_of_day"]))
    tenth_s = int(t.microsecond/100000)
    enddate = '%d%02d%02d' % (t.year,t.month,t.day)
    endtime = '%02d%02d%02d%01d' % (t.hour,t.minute,t.second,tenth_s)

    gac_io.avhrrGAC_io(satellite_name, startdate, enddate, starttime, endtime, arrLat_full, arrLon_full, ref1, ref2, ref3, bt3, bt4, bt5, arrSZA, arrSTZ, arrSAA, arrSTA, arrRAA)


if __name__ == "__main__":
    main(filename)












