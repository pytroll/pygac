#!/usr/bin/python
# Copyright (c) 2014, 2015.
#

# Author(s):
#   Sajid Pareeth <sajid.pareeth@fmach.it>
#   Martin Raspaud <martin.raspaud@smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

"""
"""

import logging
import datetime

import numpy as np

from pygac.reader import LACReader
import pygac.gac_lac_geotiepoints as gtp
import gac_io

LOG = logging.getLogger(__name__)

# LAC header object
header = np.dtype([("data_set_creation_site_id", "S3"),
                   ("ascii_blank_=_x20", "S1"),
                   ("noaa_level_1b_format_version_number", ">u2"),
                   ("noaa_level_1b_format_version_year", ">u2"),
                   ("noaa_level_1b_format_version_day_of_year", ">u2"),
                   ("reserved_for_logical_record_length", ">u2"),
                   ("reserved_for_block_size", ">u2"),
                   ("count_of_header_records", ">u2"),
                   ("zero_fill0", ">i2", (3,)),
                   ("data_set_name", "S42"),
                   ("processing_block_identification", "S8"),
                   ("noaa_spacecraft_identification_code", ">u2"),
                   ("instrument_id", ">u2"),
                   ("data_type_code", ">u2"),
                   ("tip_source_code", ">u2"),
                   ("start_of_data_set_day_count_starting_from_0_at_00h,_1_jan_1950",
                    ">u4"),
                   ("start_of_data_set_year", ">u2"),
                   ("start_of_data_set_day_of_year", ">u2"),
                   ("start_of_data_set_utc_time_of_day", ">u4"),
                   ("end_of_data_set_day_count_starting_from_0_at_00h,_1_jan_1950",
                    ">u4"),
                   ("end_of_data_set_year", ">u2"),
                   ("end_of_data_set_day_of_year", ">u2"),
                   ("end_of_data_set_utc_time_of_day", ">u4"),
                   ("year_of_last_cpids_update", ">u2"),
                   ("day_of_year_of_last_cpids_update", ">u2"),
                   ("zero_fill1", ">i2", (4,)),
                   # DATA SET QUALITY INDICATORS
                   ("instrument_status", ">u4"),
                   ("zero_fill2", ">i2"),
                   ("record_number_of_status_change", ">u2"),
                   ("second_instrument_status", ">u4"),
                   ("count_of_data_records", ">u2"),
                   ("count_of_calibrated,_earth_located_scan_lines", ">u2"),
                   ("count_of_missing_scan_lines", ">u2"),
                   ("count_of_data_gaps", ">u2"),
                   ("count_of_data_frames_without_frame_sync_word", ">u2"),
                   ("count_of_pacs_detected_tip_parity_errors", ">u2"),
                   ("sum_of_all_auxiliary_sync_errors_detected", ">u2"),
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
                   ("zero_fill4", ">i2", (5,)),
                   # CALIBRATION
                   ("ramp_auto_calibration_indicators_bit_field", ">u2"),
                   ("year_of_most_recent_solar_channel_calibration", ">u2"),
                   ("day_of_year_of_most_recent_solar_channel_calibration",
                    ">u2"),
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
                   ("zero_fill5", ">i4", (2,)),
                   # RADIANCE CONVERSION
                   ("ch_1_solar_filtered_irradiance_in_wavelength", ">i4"),
                   ("ch_1_equivalent_filter_width_in_wavelength", ">i4"),
                   ("ch_2_solar_filtered_irradiance_in_wavelength", ">i4"),
                   ("ch_2_equivalent_filter_width_in_wavelength", ">i4"),
                   ("ch_3a_solar_filtered_irradiance_in_wavelength", ">i4"),
                   ("ch_3a_equivalent_filter_width_in_wavelength", ">i4"),
                   ("ch_3b_central_wavenumber", ">i4"),
                   ("ch_3b_constant_1", ">i4"),
                   ("ch_3b_constant_2", ">i4"),
                   ("ch_4_central_wavenumber", ">i4"),
                   ("ch_4_constant_1", ">i4"),
                   ("ch_4_constant_2", ">i4"),
                   ("ch_5_central_wavenumber", ">i4"),
                   ("ch_5_constant_1", ">i4"),
                   ("ch_5_constant_2", ">i4"),
                   ("zero_fill6", ">i4", (3,)),
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
                   ("zero_fill8", ">i4", (4,)),
                   # ANALOG TELEMETRY CONVERSION
                   ("patch_temperature_conversion_coefficient_1", ">i4"),
                   ("patch_temperature_conversion_coefficient_2", ">i4"),
                   ("patch_temperature_conversion_coefficient_3", ">i4"),
                   ("patch_temperature_conversion_coefficient_4", ">i4"),
                   ("patch_temperature_conversion_coefficient_5", ">i4"),
                   ("patch_temperature_conversion_coefficient_6", ">i4"),
                   ("patch_temperature_extended_conversion_coefficient_1",
                    ">i4"),
                   ("patch_temperature_extended_conversion_coefficient_2",
                    ">i4"),
                   ("patch_temperature_extended_conversion_coefficient_3",
                    ">i4"),
                   ("patch_temperature_extended_conversion_coefficient_4",
                    ">i4"),
                   ("patch_temperature_extended_conversion_coefficient_5",
                    ">i4"),
                   ("patch_temperature_extended_conversion_coefficient_6", ">i4"),
                   ("patch_power_conversion_coefficient_1", ">i4"),
                   ("patch_power_conversion_coefficient_2", ">i4"),
                   ("patch_power_conversion_coefficient_3", ">i4"),
                   ("patch_power_conversion_coefficient_4", ">i4"),
                   ("patch_power_conversion_coefficient_5", ">i4"),
                   ("patch_power_conversion_coefficient_6", ">i4"),
                   ("radiator_temperature_conversion_coefficient_1", ">i4"),
                   ("radiator_temperature_conversion_coefficient_2", ">i4"),
                   ("radiator_temperature_conversion_coefficient_3", ">i4"),
                   ("radiator_temperature_conversion_coefficient_4", ">i4"),
                   ("radiator_temperature_conversion_coefficient_5", ">i4"),
                   ("radiator_temperature_conversion_coefficient_6", ">i4"),
                   ("blackbody_temperature_1_conversion_coefficient_1", ">i4"),
                   ("blackbody_temperature_1_conversion_coefficient_2", ">i4"),
                   ("blackbody_temperature_1_conversion_coefficient_3", ">i4"),
                   ("blackbody_temperature_1_conversion_coefficient_4", ">i4"),
                   ("blackbody_temperature_1_conversion_coefficient_5", ">i4"),
                   ("blackbody_temperature_1_conversion_coefficient_6", ">i4"),
                   ("blackbody_temperature_2_conversion_coefficient_1", ">i4"),
                   ("blackbody_temperature_2_conversion_coefficient_2", ">i4"),
                   ("blackbody_temperature_2_conversion_coefficient_3", ">i4"),
                   ("blackbody_temperature_2_conversion_coefficient_4", ">i4"),
                   ("blackbody_temperature_2_conversion_coefficient_5", ">i4"),
                   ("blackbody_temperature_2_conversion_coefficient_6", ">i4"),
                   ("blackbody_temperature_3_conversion_coefficient_1", ">i4"),
                   ("blackbody_temperature_3_conversion_coefficient_2", ">i4"),
                   ("blackbody_temperature_3_conversion_coefficient_3", ">i4"),
                   ("blackbody_temperature_3_conversion_coefficient_4", ">i4"),
                   ("blackbody_temperature_3_conversion_coefficient_5", ">i4"),
                   ("blackbody_temperature_3_conversion_coefficient_6", ">i4"),
                   ("blackbody_temperature_4_conversion_coefficient_1", ">i4"),
                   ("blackbody_temperature_4_conversion_coefficient_2", ">i4"),
                   ("blackbody_temperature_4_conversion_coefficient_3", ">i4"),
                   ("blackbody_temperature_4_conversion_coefficient_4", ">i4"),
                   ("blackbody_temperature_4_conversion_coefficient_5", ">i4"),
                   ("blackbody_temperature_4_conversion_coefficient_6", ">i4"),
                   ("electronics_current_conversion_coefficient_1", ">i4"),
                   ("electronics_current_conversion_coefficient_2", ">i4"),
                   ("electronics_current_conversion_coefficient_3", ">i4"),
                   ("electronics_current_conversion_coefficient_4", ">i4"),
                   ("electronics_current_conversion_coefficient_5", ">i4"),
                   ("electronics_current_conversion_coefficient_6", ">i4"),
                   ("motor_current_conversion_coefficient_1", ">i4"),
                   ("motor_current_conversion_coefficient_2", ">i4"),
                   ("motor_current_conversion_coefficient_3", ">i4"),
                   ("motor_current_conversion_coefficient_4", ">i4"),
                   ("motor_current_conversion_coefficient_5", ">i4"),
                   ("motor_current_conversion_coefficient_6", ">i4"),
                   ("earth_shield_position_conversion_coefficient_1", ">i4"),
                   ("earth_shield_position_conversion_coefficient_2", ">i4"),
                   ("earth_shield_position_conversion_coefficient_3", ">i4"),
                   ("earth_shield_position_conversion_coefficient_4", ">i4"),
                   ("earth_shield_position_conversion_coefficient_5", ">i4"),
                   ("earth_shield_position_conversion_coefficient_6", ">i4"),
                   ("electronics_temperature_conversion_coefficient_1", ">i4"),
                   ("electronics_temperature_conversion_coefficient_2", ">i4"),
                   ("electronics_temperature_conversion_coefficient_3", ">i4"),
                   ("electronics_temperature_conversion_coefficient_4", ">i4"),
                   ("electronics_temperature_conversion_coefficient_5", ">i4"),
                   ("electronics_temperature_conversion_coefficient_6", ">i4"),
                   ("cooler_housing_temperature_conversion_coefficient_1", ">i4"),
                   ("cooler_housing_temperature_conversion_coefficient_2", ">i4"),
                   ("cooler_housing_temperature_conversion_coefficient_3", ">i4"),
                   ("cooler_housing_temperature_conversion_coefficient_4", ">i4"),
                   ("cooler_housing_temperature_conversion_coefficient_5", ">i4"),
                   ("cooler_housing_temperature_conversion_coefficient_6", ">i4"),
                   ("baseplate_temperature_conversion_coefficient_1", ">i4"),
                   ("baseplate_temperature_conversion_coefficient_2", ">i4"),
                   ("baseplate_temperature_conversion_coefficient_3", ">i4"),
                   ("baseplate_temperature_conversion_coefficient_4", ">i4"),
                   ("baseplate_temperature_conversion_coefficient_5", ">i4"),
                   ("baseplate_temperature_conversion_coefficient_6", ">i4"),
                   ("motor_housing_temperature_conversion_coefficient_1", ">i4"),
                   ("motor_housing_temperature_conversion_coefficient_2", ">i4"),
                   ("motor_housing_temperature_conversion_coefficient_3", ">i4"),
                   ("motor_housing_temperature_conversion_coefficient_4", ">i4"),
                   ("motor_housing_temperature_conversion_coefficient_5", ">i4"),
                   ("motor_housing_temperature_conversion_coefficient_6", ">i4"),
                   ("a/d_converter_temperature_conversion_coefficient_1", ">i4"),
                   ("a/d_converter_temperature_conversion_coefficient_2", ">i4"),
                   ("a/d_converter_temperature_conversion_coefficient_3", ">i4"),
                   ("a/d_converter_temperature_conversion_coefficient_4", ">i4"),
                   ("a/d_converter_temperature_conversion_coefficient_5", ">i4"),
                   ("a/d_converter_temperature_conversion_coefficient_6", ">i4"),
                   ("detector_#4_bias_voltage_conversion_coefficient_1", ">i4"),
                   ("detector_#4_bias_voltage_conversion_coefficient_2", ">i4"),
                   ("detector_#4_bias_voltage_conversion_coefficient_3", ">i4"),
                   ("detector_#4_bias_voltage_conversion_coefficient_4", ">i4"),
                   ("detector_#4_bias_voltage_conversion_coefficient_5", ">i4"),
                   ("detector_#4_bias_voltage_conversion_coefficient_6", ">i4"),
                   ("detector_#5_bias_voltage_conversion_coefficient_1", ">i4"),
                   ("detector_#5_bias_voltage_conversion_coefficient_2", ">i4"),
                   ("detector_#5_bias_voltage_conversion_coefficient_3", ">i4"),
                   ("detector_#5_bias_voltage_conversion_coefficient_4", ">i4"),
                   ("detector_#5_bias_voltage_conversion_coefficient_5", ">i4"),
                   ("detector_#5_bias_voltage_conversion_coefficient_6", ">i4"),
                   ("blackbody_temperature_channel3b_conversion_coefficient_1", ">i4"),
                   ("blackbody_temperature_channel3b_conversion_coefficient_2", ">i4"),
                   ("blackbody_temperature_channel3b_conversion_coefficient_3", ">i4"),
                   ("blackbody_temperature_channel3b_conversion_coefficient_4", ">i4"),
                   ("blackbody_temperature_channel3b_conversion_coefficient_5", ">i4"),
                   ("blackbody_temperature_channel3b_conversion_coefficient_6", ">i4"),
                   ("blackbody_temperature_channel4_conversion_coefficient_1",
                    ">i4"),
                   ("blackbody_temperature_channel4_conversion_coefficient_2",
                    ">i4"),
                   ("blackbody_temperature_channel4_conversion_coefficient_3",
                    ">i4"),
                   ("blackbody_temperature_channel4_conversion_coefficient_4",
                    ">i4"),
                   ("blackbody_temperature_channel4_conversion_coefficient_5",
                    ">i4"),
                   ("blackbody_temperature_channel4_conversion_coefficient_6",
                    ">i4"),
                   ("blackbody_temperature_channel5_conversion_coefficient_1",
                    ">i4"),
                   ("blackbody_temperature_channel5_conversion_coefficient_2",
                    ">i4"),
                   ("blackbody_temperature_channel5_conversion_coefficient_3",
                    ">i4"),
                   ("blackbody_temperature_channel5_conversion_coefficient_4",
                    ">i4"),
                   ("blackbody_temperature_channel5_conversion_coefficient_5",
                    ">i4"),
                   ("blackbody_temperature_channel5_conversion_coefficient_6",
                    ">i4"),
                   ("reference_voltage_conversion_coefficient_1", ">i4"),
                   ("reference_voltage_conversion_coefficient_2", ">i4"),
                   ("reference_voltage_conversion_coefficient_3", ">i4"),
                   ("reference_voltage_conversion_coefficient_4", ">i4"),
                   ("reference_voltage_conversion_coefficient_5", ">i4"),
                   ("reference_voltage_conversion_coefficient_6", ">i4"),
                   # METOP MANEUVERS IDENTIFICATION
                   ("start_of_maneuver_year", ">u2"),
                   ("start_of_maneuver_day", ">u2"),
                   ("start_of_maneuver_utc_time_of_day", ">u4"),
                   ("end_of_manuever_year", ">u2"),
                   ("end_of_maneuver_day", ">u2"),
                   ("end_of_maneuver_utc_time_of_day", ">u4"),
                   ("change_in_spacecraft_velocity", ">i4"),
                   ("change in spacecraft_mass", ">u4"),
                   # CLOUDS FROM AVHRR(CLAVR)
                   ("clavr_status_bit_field", ">u2"),
                   ("zero_fill9", ">i2")])
# FILLER
#                  ("zero_fill10", ">i4", (4608 - 688, ))])


# video data object
scanline = np.dtype([("scan_line_number", ">u2"),
                     ("scan_line_year", ">i2"),
                     ("scan_line_day_of_year", ">u2"),
                     ("satellite_clock_drift_delta", ">i2"),
                     ("scan_line_utc_time_of_day", ">u4"),
                     ("scan_line_bit_field", ">u2"),
                     ("zero_fill0", ">i2", (5,)),
                     # QUALITY INDICATORS
                     ("quality_indicator_bit_field", ">u4"),
                     ("scan_line_quality_flags_reserved", ">u1"),
                     ("scan_line_quality_flags_time_problem_code", ">u1"),
                     ("scan_line_quality_flags_calibration_problem_code", ">u1"),
                     ("scan_line_quality_flags_earth_location_problem_code", ">u1"),
                     ("calibration_quality_flags", ">u2", (3,)),
                     ("count_of_bit_errors_in_frame_sync", ">u2"),
                     ("zero_fill1", ">i4", (2,)),
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
                     ("computed_yaw_steering", ">i2", (3,)),
                     ("total_applied_attitude_correction", ">i2", (3,)),
                     ("navigation_status_bit_field", ">u4"),
                     ("time_associated_with_euler_angles", ">u4"),
                     ("euler_angles", ">i2", (3,)),
                     ("spacecraft_altitude_above_reference_ellipsoid", ">u2"),
                     ("angular_relationships", ">i2", (153,)),
                     ("zero_fill2", ">i2", (3,)),
                     ("earth_location", ">i4", (102,)),
                     ("zero_fill3", ">i4", (2,)),
                     # FRAME TELEMETRY
                     ("frame_sync", ">u2", (6,)),
                     ("id", ">u2", (2,)),
                     ("time_code", ">u2", (4,)),
                     ("ramp_calibration", ">u2", (5,)),
                     ("internal_target_temperature", ">u2", (3,)),
                     ("patch_temperature", ">u2"),
                     ("undefined", ">u2"),
                     ("back_scan", ">u2", (30,)),
                     ("space_data", ">u2", (50,)),
                     ("sync_delta", ">u2"),
                     ("zero_fill4", ">i2"),
                     # EARTH OBSERVATIONS
                     ("sensor_data", ">u4", (3414,)),
                     ("zero_fill5", ">i4", (2,)),
                     # DIGITAL B HOUSEKEEPING TELEMETRY
                     ("digital_b_telemetry_update_flags", ">u2"),
                     ("avhrr_digital_b_data", ">u2"),
                     ("zero_fill6", ">i4", (3,)),
                     # ANALOG HOUSEKEEPING DATA (TIP)
                     ("analog_telemetry_update_flags", ">u4"),
                     ("patch_temperature_range_0_255", ">u1"),
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
                     ("zero_fill7", ">i2", (3,)),
                     # CLOUDS FROM AVHRR (CLAVR)
                     ("reserved_clavr_status_bit_field", ">u4"),
                     ("reserved_clavr", ">u4"),
                     ("reserved_clavr_ccm", ">u2", (256,)),
                     # FILLER
                     ("zero_fill8", ">i4", (94,))])


class KLMReader(LACReader):
    instrument_ids = {4:  15,
                      2:  16,
                      6:  17,
                      7:  18,
                      8:  19,
                      12: 'a',
                      11: 'b',
                      }

    spacecraft_names = {4:  'noaa15',
                        2:  'noaa16',
                        6:  'noaa17',
                        7:  'noaa18',
                        8:  'noaa19',
                        12: 'metopa',
                        11: 'metopb',
                        }
    spacecrafts_orbital = {4:  'noaa 15',
                           2:  'noaa 16',
                           6:  'noaa 17',
                           7:  'noaa 18',
                           8:  'noaa 19',
                           12: 'metop 02',
                           11: 'metop 01',
                           }

    def read(self, filename):
        with open(filename) as fd_:
            self.head = np.fromfile(fd_, dtype=header, count=1)[0]
            # The value below is very important, in this case:15872, from the LAC header table last row.
            fd_.seek(15872, 0)
            self.scans = np.fromfile(
                fd_, dtype=scanline, count=self.head["count_of_data_records"])

        self.spacecraft_id = self.head["noaa_spacecraft_identification_code"]
        self.instrument_id = self.instrument_ids[self.spacecraft_id]
        self.spacecraft_name = self.spacecraft_names[self.spacecraft_id]

    def get_telemetry(self):
        prt_counts = np.mean(self.scans["internal_target_temperature"][:, 0:3], axis=1)

        # getting ICT counts

        ict_counts = np.zeros((len(self.scans), 3))
        ict_counts[:, 0] = np.mean(self.scans["back_scan"][:, 0::3], axis=1)
        ict_counts[:, 1] = np.mean(self.scans["back_scan"][:, 1::3], axis=1)
        ict_counts[:, 2] = np.mean(self.scans["back_scan"][:, 2::3], axis=1)

        # getting space counts

        space_counts = np.zeros((len(self.scans), 3))
        space_counts[:, 0] = np.mean(self.scans["space_data"][:, 2::5], axis=1)
        space_counts[:, 1] = np.mean(self.scans["space_data"][:, 3::5], axis=1)
        space_counts[:, 2] = np.mean(self.scans["space_data"][:, 4::5], axis=1)

        return prt_counts, ict_counts, space_counts

    def get_lonlat(self):
        arr_lat = self.scans["earth_location"][:, 0::2] / 1e4
        arr_lon = self.scans["earth_location"][:, 1::2] / 1e4

        self.lons, self.lats = gtp.lac_geo_interpolator(arr_lon, arr_lat)
        return self.lons, self.lats

    def get_times(self):
        if self.utcs is None:
            year = self.scans["scan_line_year"]
            jday = self.scans["scan_line_day_of_year"]
            msec = self.scans["scan_line_utc_time_of_day"]

            self.utcs = (((year - 1970).astype('datetime64[Y]')
                          + (jday - 1).astype('timedelta64[D]')).astype('datetime64[ms]')
                         + msec.astype('timedelta64[ms]'))
            self.times = self.utcs.astype(datetime.datetime)
        return self.utcs

    def get_ch3_switch(self):
        """0 for 3b, 1 for 3a
        """
        return self.scans["scan_line_bit_field"][:] & 3

    def get_corrupt_mask(self):
        # corrupt scanlines

        mask = ((self.scans["quality_indicator_bit_field"] >> 31) |
                ((self.scans["quality_indicator_bit_field"] << 3) >> 31) |
                ((self.scans["quality_indicator_bit_field"] << 4) >> 31))

        number_of_scans = self.scans["internal_target_temperature"].shape[0]
        qual_flags = np.zeros((int(number_of_scans), 7))
        qual_flags[:, 0] = self.scans["scan_line_number"]
        qual_flags[:, 1] = (self.scans["quality_indicator_bit_field"] >> 31)
        qual_flags[:, 2] = ((self.scans["quality_indicator_bit_field"] << 3) >> 31)
        qual_flags[:, 3] = ((self.scans["quality_indicator_bit_field"] << 4) >> 31)
        qual_flags[:, 4] = ((self.scans["quality_indicator_bit_field"] << 24) >> 30)
        qual_flags[:, 5] = ((self.scans["quality_indicator_bit_field"] << 26) >> 30)
        qual_flags[:, 6] = ((self.scans["quality_indicator_bit_field"] << 28) >> 30)

        return mask.astype(bool), qual_flags


def main(filename):
    tic = datetime.datetime.now()
    reader = KLMReader()
    reader.read(filename)
    reader.get_lonlat()
    reader.adjust_clock_drift()
    channels = reader.get_calibrated_channels()
    sat_azi, sat_zen, sun_azi, sun_zen, rel_azi = reader.get_angles()

    mask = reader.get_corrupt_mask()
    gac_io.save_gac(reader.spacecraft_name,
                    reader.times[0], reader.times[1],
                    reader.lats, reader.lons,
                    channels[:, :, 0], channels[:, :, 1],
                    np.ones_like(channels[:, :, 0]) * -1,
                    channels[:, :, 2],
                    channels[:, :, 3],
                    channels[:, :, 4],
                    sun_zen, sat_zen, sun_azi, sat_azi, rel_azi,
                    mask)
    LOG.info("pygac took: %s", str(datetime.datetime.now() - tic))


if __name__ == "__main__":
    import sys

    main(sys.argv[1])
