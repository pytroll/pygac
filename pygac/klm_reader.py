#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014-2019

# Author(s):

#   Abhay Devasthale <abhay.devasthale@smhi.se>
#   Sajid Pareeth <sajid.pareeth@fmach.it>
#   Martin Raspaud <martin.raspaud@smhi.se>
#   Adam Dybbroe <adam.dybbroe@smhi.se>
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

"""Read KLM data.

Reads L1b GAC/LAC data from KLM series of satellites (NOAA-15 and later).
Format specification can be found in section 8 of the `KLM user guide`_.

.. _KLM user guide:
    https://www.ncei.noaa.gov/pub/data/satellite/publications/podguides/N-15%20thru%20N-19/

"""

import datetime
import logging
try:
    from enum import IntFlag
except ImportError:
    # python version < 3.6, use a simple object without nice representation
    IntFlag = object

import numpy as np

from pygac.correct_tsm_issue import TSM_AFFECTED_INTERVALS_KLM, get_tsm_idx
from pygac.reader import Reader, ReaderError
from pygac.utils import file_opener

LOG = logging.getLogger(__name__)


class KLM_QualityIndicator(IntFlag):
    """Quality Indicators.

    Source: KLM guide

    - Table 8.3.1.3.3.1-1. Format of packed LAC/HRPT Data Sets (Version 2,
      pre-April 28, 2005).
    - Table 8.3.1.3.3.2-1. Format of LAC/HRPT Data Record for NOAA-N
      (Version 5, post-November 14, 2006, all spacecraft).
    - Table 8.3.1.4.3.1-1. Format of packed GAC Data Record for NOAA KLM
      (Version 2, pre-April 28, 2005).
    - Table 8.3.1.4.3.2-1. Format of GAC Data Record for NOAA-N (Version 4,
      post-January 25, 2006, all spacecraft).

    Notes:
    - Table 8.3.1.3.3.1-1. and Table 8.3.1.4.3.1-1. define bit: 21 as
      "frame sync word not valid"
    - Table 8.3.1.3.3.2-1. and Table 8.3.1.4.3.2-1. define bit: 21 as
      "flywheeling detected during this frame"

    """

    FATAL_FLAG = 2**31  # Data should not be used for product generation
    TIME_ERROR = 2**30  # Time sequence error detected within this scan
    DATA_GAP = 2**29  # Data gap precedes this scan
    CALIBRATION = 2**28  # Insufficient data for calibration
    NO_EARTH_LOCATION = 2**27  # Earth location data not available
    CLOCK_UPDATE = 2**26  # First good time following a clock update (nominally 0)
    INSTRUMENT_CHANGE = 2**25  # Instrument status changed with this scan
    BIT_SYNC_STATUS = 2**24  # Sync lock dropped during this frame
    SYNC_ERROR = 2**23  # Frame sync word error greater than zero
    FRAME_SYNC_LOCK = 2**22  # Frame sync previously dropped lock
    SYNC_INVALID = 2**21  # Frame sync word not valid
    FLYWHEELING = 2**21  # Flywheeling detected during this frame
    BIT_SLIPPAGE = 2**20  # Bit slippage detected during this frame
    # Note: Bit 19 - 9 are not defined for KLMs
    TIP_PARITY = 2**8  # TIP parity error detected
    # Reflected Sunlight (RS) detected (solar blackbody contamination)
    CH_3B_RS = 2**7
    CH_3B_RS_ANOMALY = 2**6
    CH_3_CONTAMINATION = CH_3B_RS | CH_3B_RS_ANOMALY  # POD compatible alias
    CH_4_RS = 2**5
    CH_4_RS_ANOMALY = 2**4
    CH_4_CONTAMINATION = CH_4_RS | CH_4_RS_ANOMALY  # POD compatible alias
    CH_5_RS = 2**3
    CH_5_RS_ANOMALY = 2**2
    CH_5_CONTAMINATION = CH_5_RS | CH_5_RS_ANOMALY  # POD compatible alias
    DATA_JITTER = 2**1  # Resync occurred on this frame
    PSEUDO_NOISE = 2**0  # Pseudo noise occurred on this frame


# header object

header = np.dtype([("data_set_creation_site_id", "S3"),
                   ("ascii_blank_=_x20", "S1"),
                   ("noaa_level_1b_format_version_number", ">u2"),
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
                   # only in version 5
                   ("count_of_data_frames_without_frame_sync_word", ">u2"),
                   ("count_of_pacs_detected_tip_parity_errors", ">u2"),
                   # only in version 5
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
                   ("zero_fill4", ">i2", (5, )),
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
                   ("zero_fill5", ">i4", (2, )),
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
                   ("zero_fill8", ">i4", (4, ))])

# ANALOG TELEMETRY CONVERSION
analog_telemetry_v2 = np.dtype([("patch_temperature_conversion_coefficient_1", ">i2"),
                                ("patch_temperature_conversion_coefficient_2", ">i2"),
                                ("patch_temperature_conversion_coefficient_3", ">i2"),
                                ("patch_temperature_conversion_coefficient_4", ">i2"),
                                ("patch_temperature_conversion_coefficient_5", ">i2"),
                                ("reserved0", ">i2"),
                                ("patch_temperature_extended_conversion_coefficient_1",
                                 ">i2"),
                                ("patch_temperature_extended_conversion_coefficient_2",
                                 ">i2"),
                                ("patch_temperature_extended_conversion_coefficient_3",
                                 ">i2"),
                                ("patch_temperature_extended_conversion_coefficient_4",
                                 ">i2"),
                                ("patch_temperature_extended_conversion_coefficient_5",
                                 ">i2"),
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
                                ("radiato   r_temperature_conversion_coefficient_5", ">i2"),
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
                                ("cooler_housing_temperature_conversion_coefficient_1",
                                 ">i2"),
                                ("cooler_housing_temperature_conversion_coefficient_2",
                                 ">i2"),
                                ("cooler_housing_temperature_conversion_coefficient_3",
                                 ">i2"),
                                ("cooler_housing_temperature_conversion_coefficient_4",
                                 ">i2"),
                                ("cooler_housing_temperature_conversion_coefficient_5",
                                 ">i2"),
                                ("reserved12", ">i2"),
                                ("baseplate_temperature_conversion_coefficient_1", ">i2"),
                                ("baseplate_temperature_conversion_coefficient_2", ">i2"),
                                ("baseplate_temperature_conversion_coefficient_3", ">i2"),
                                ("baseplate_temperature_conversion_coefficient_4", ">i2"),
                                ("baseplate_temperature_conversion_coefficient_5", ">i2"),
                                ("reserved13", ">i2"),
                                ("motor_housing_temperature_conversion_coefficient_1",
                                 ">i2"),
                                ("motor_housing_temperature_conversion_coefficient_2",
                                 ">i2"),
                                ("motor_housing_temperature_conversion_coefficient_3",
                                 ">i2"),
                                ("motor_housing_temperature_conversion_coefficient_4",
                                 ">i2"),
                                ("motor_housing_temperature_conversion_coefficient_5",
                                 ">i2"),
                                ("reserved14", ">i2"),
                                ("a/d_converter_temperature_conversion_coefficient_1",
                                 ">i2"),
                                ("a/d_converter_temperature_conversion_coefficient_2",
                                 ">i2"),
                                ("a/d_converter_temperature_conversion_coefficient_3",
                                 ">i2"),
                                ("a/d_converter_temperature_conversion_coefficient_4",
                                 ">i2"),
                                ("a/d_converter_temperature_conversion_coefficient_5",
                                 ">i2"),
                                ("reserved15", ">i2"),
                                ("detector_#4_bias_voltage_conversion_coefficient_1",
                                 ">i2"),
                                ("detector_#4_bias_voltage_conversion_coefficient_2",
                                 ">i2"),
                                ("detector_#4_bias_voltage_conversion_coefficient_3",
                                 ">i2"),
                                ("detector_#4_bias_voltage_conversion_coefficient_4",
                                 ">i2"),
                                ("detector_#4_bias_voltage_conversion_coefficient_5",
                                 ">i2"),
                                ("reserved16", ">i2"),
                                ("detector_#5_bias_voltage_conversion_coefficient_1",
                                 ">i2"),
                                ("detector_#5_bias_voltage_conversion_coefficient_2",
                                 ">i2"),
                                ("detector_#5_bias_voltage_conversion_coefficient_3",
                                 ">i2"),
                                ("detector_#5_bias_voltage_conversion_coefficient_4",
                                 ">i2"),
                                ("detector_#5_bias_voltage_conversion_coefficient_5",
                                 ">i2"),
                                ("reserved17", ">i2"),
                                ("channel_3b_blackbody_view_conversion_coefficient_1",
                                 ">i2"),
                                ("channel_3b_blackbody_view_conversion_coefficient_2",
                                 ">i2"),
                                ("channel_3b_blackbody_view_conversion_coefficient_3",
                                 ">i2"),
                                ("channel_3b_blackbody_view_conversion_coefficient_4",
                                 ">i2"),
                                ("channel_3b_blackbody_view_conversion_coefficient_5",
                                 ">i2"),
                                ("reserved18", ">i2"),
                                ("channel_4_blackbody_view_conversion_coefficient_1",
                                 ">i2"),
                                ("channel_4_blackbody_view_conversion_coefficient_2",
                                 ">i2"),
                                ("channel_4_blackbody_view_conversion_coefficient_3",
                                 ">i2"),
                                ("channel_4_blackbody_view_conversion_coefficient_4",
                                 ">i2"),
                                ("channel_4_blackbody_view_conversion_coefficient_5",
                                 ">i2"),
                                ("reserved19", ">i2"),
                                ("channel_5_blackbody_view_conversion_coefficient_1",
                                 ">i2"),
                                ("channel_5_blackbody_view_conversion_coefficient_2",
                                 ">i2"),
                                ("channel_5_blackbody_view_conversion_coefficient_3",
                                 ">i2"),
                                ("channel_5_blackbody_view_conversion_coefficient_4",
                                 ">i2"),
                                ("channel_5_blackbody_view_conversion_coefficient_5",
                                 ">i2"),
                                ("reserved20", ">i2"),
                                ("reference_voltage_conversion_coefficient_1", ">i2"),
                                ("reference_voltage_conversion_coefficient_2", ">i2"),
                                ("reference_voltage_conversion_coefficient_3", ">i2"),
                                ("reference_voltage_conversion_coefficient_4", ">i2"),
                                ("reference_voltage_conversion_coefficient_5", ">i2"),
                                ("reserved21", ">i2")])
# FILLER
#                  ("zero_fill9", ">i4", (4608 - 688, ))])


analog_telemetry_v5 = np.dtype([("patch_temperature_conversion_coefficient_1", ">i4"),
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


ars_header = np.dtype([('COST_number', 'S6'),
                       ('SAA_number', 'S8'),
                       ('order_creation_year', 'S4'),
                       ('order_creation_day_of_year', 'S3'),
                       ('processing_site_code', 'S1'),
                       ('processing_software', 'S8'),
                       # data selection criteria
                       ('data_set_name', 'S42'),
                       ('ascii_blank_', 'S2'),
                       ('select_flag', 'S1'),
                       ('beginning_latitude', 'S3'),
                       ('ending_latitude', 'S3'),
                       ('beginning_longitude', 'S4'),
                       ('ending_longitude', 'S4'),
                       ('start_hour', 'S2'),
                       ('start_minute', 'S2'),
                       ('number_of_minutes', 'S3'),
                       ('appended_data_flag', 'S1'),
                       ('channel_select_flag', 'S1', (20, )),
                       # dataset summary
                       ('ascii_blank__', 'S29'),
                       ('ascend_descend_flag', 'S1'),
                       ('first_latitude', 'S3'),
                       ('last_latitude', 'S3'),
                       ('first_longitude', 'S4'),
                       ('last_longitude', 'S4'),
                       ('data_format', 'S20'),
                       ('size_of_record', 'S6'),
                       ('number_of_records', 'S6'),
                       # filler
                       ('ascii_blank', 'S319')
                       ])


class KLMReader(Reader):
    """Reader for KLM data."""

    spacecraft_names = {4: 'noaa15',
                        2: 'noaa16',
                        6: 'noaa17',
                        7: 'noaa18',
                        8: 'noaa19',
                        12: 'metopa',
                        11: 'metopb',
                        13: 'metopc',
                        }
    spacecrafts_orbital = {4: 'noaa 15',
                           2: 'noaa 16',
                           6: 'noaa 17',
                           7: 'noaa 18',
                           8: 'noaa 19',
                           12: 'metop 02',
                           11: 'metop 01',
                           13: 'metop 03',
                           }

    tsm_affected_intervals = TSM_AFFECTED_INTERVALS_KLM

    QFlag = KLM_QualityIndicator
    _quality_indicators_key = "quality_indicator_bit_field"

    def read(self, filename, fileobj=None):
        """Read the data.

        Args:
            filename: Path to GAC/LAC file
            fileobj: An open file object to read from. (optional)

        Returns:
            header: numpy record array
                The header metadata
            scans: numpy record array
                The scanlines

        """
        # Note that np.fromfile does not work with gzip.GzipFile
        # objects (numpy version 1.16.4), because it restricts the
        # file objects to (io.FileIO, io.BufferedReader, io.BufferedWriter)
        # see: numpy.compat.py3k.isfileobj
        self.filename = filename
        LOG.info('Reading %s', self.filename)
        with file_opener(fileobj or filename) as fd_:
            self.ars_head, self.head = self.read_header(
                filename, fileobj=fd_)
            if self.ars_head:
                ars_offset = ars_header.itemsize
            else:
                ars_offset = 0
            self.header_version = self.head[
                "noaa_level_1b_format_version_number"]
            if self.header_version >= 5:
                self.analog_telemetry, = np.frombuffer(
                    fd_.read(analog_telemetry_v5.itemsize),
                    dtype=analog_telemetry_v5, count=1)
            else:
                self.analog_telemetry, = np.frombuffer(
                    fd_.read(analog_telemetry_v2.itemsize),
                    dtype=analog_telemetry_v2, count=1)
            # LAC: 1, GAC: 2, ...
            self.data_type = self.head['data_type_code']
            # read until end of file
            fd_.seek(self.offset + ars_offset, 0)
            buffer = fd_.read()
            count = self.head["count_of_data_records"]
            self._read_scanlines(buffer, count)
        self.correct_scan_line_numbers()
        self.spacecraft_id = self.head["noaa_spacecraft_identification_code"]
        self.spacecraft_name = self.spacecraft_names[self.spacecraft_id]

    @classmethod
    def read_header(cls, filename, fileobj=None):
        """Read the file header.

        Args:
            filename (str): Path to GAC/LAC file
            fileobj: An open file object to read from. (optional)

        Returns:
            archive_header (struct): archive header
            header (struct): file header
        """
        with file_opener(fileobj or filename) as fd_:
            # read ars_header if present
            _ars_head, = np.frombuffer(
                fd_.read(ars_header.itemsize),
                dtype=ars_header, count=1)
            if _ars_head['data_format'].startswith(b'NOAA Level 1b'):
                ars_head = _ars_head.copy()
            else:
                fd_.seek(0)
                ars_head = None
            # need to copy frombuffer to have write access on head
            head, = np.frombuffer(
                fd_.read(header.itemsize),
                dtype=header, count=1).copy()
        head = cls._correct_data_set_name(head, filename)
        cls._validate_header(head)
        return ars_head, head

    @classmethod
    def _validate_header(cls, header):
        """Check if the header belongs to this reader."""
        # call super to enter the Method Resolution Order (MRO)
        super(KLMReader, cls)._validate_header(header)
        LOG.debug("validate header")
        data_set_name = header['data_set_name'].decode()
        # split header into parts
        creation_site, transfer_mode, platform_id = (
            data_set_name.split('.')[:3])
        allowed_ids = ['NK', 'NL', 'NM', 'NN', 'NP', 'M1', 'M2', 'M3']
        if platform_id not in allowed_ids:
            raise ReaderError('Improper platform id "%s"!' % platform_id)

    def get_telemetry(self):
        """Get the telemetry.

        Returns:
            prt_counts: np.array
            ict_counts: np.array
            space_counts: np.array

        """
        prt_counts = np.mean(self.scans["telemetry"]['PRT'], axis=1)

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

    def _get_lonlat(self):
        """Get the longitudes and latitudes."""
        lats = self.scans["earth_location"][:, 0::2] / 1e4
        lons = self.scans["earth_location"][:, 1::2] / 1e4
        return lons, lats

    def get_header_timestamp(self):
        """Get the timestamp from the header.

        Returns:
            A datetime object containing the timestamp from the header.

        Raises:
            A ValueError if the timestamp is corrupt.

        """
        year = self.head['start_of_data_set_year']
        jday = self.head['start_of_data_set_day_of_year']
        msec = self.head['start_of_data_set_utc_time_of_day']
        try:
            return self.to_datetime(self.to_datetime64(year=year, jday=jday,
                                                       msec=msec))
        except ValueError as err:
            raise ValueError('Corrupt header timestamp: {0}'.format(err))

    def _get_times(self):
        """Get the times of the scanlines."""
        year = self.scans["scan_line_year"]
        jday = self.scans["scan_line_day_of_year"]
        msec = self.scans["scan_line_utc_time_of_day"]
        return year, jday, msec

    def get_ch3_switch(self):
        """Channel 3 identification.

        0: Channel 3b (Brightness temperature)
        1: Channel 3a (Reflectance)
        2: Transition (No data)
        """
        return self.scans["scan_line_bit_field"][:] & 3

    def _get_ir_channels_to_calibrate(self):
        ir_channels_to_calibrate = [3, 4, 5]
        if np.all(self.get_ch3_switch() != 0):
            ir_channels_to_calibrate = [4, 5]
        return ir_channels_to_calibrate

    def postproc(self, channels):
        """Apply KLM specific postprocessing.

        Masking out 3a/3b/transition.
        """
        switch = self.get_ch3_switch()
        channels[:, :, 2][switch == 0] = np.nan
        channels[:, :, 3][switch == 1] = np.nan
        channels[:, :, 2][switch == 2] = np.nan
        channels[:, :, 3][switch == 2] = np.nan

    def _adjust_clock_drift(self):
        """Adjust the geolocation to compensate for the clock error.

        Note:
            Clock drift correction is only applied to POD satellites.
            On the KLM series, the clock is updated daily.
        """

    def get_tsm_pixels(self, channels):
        """Determine pixels affected by the scan motor issue.

        Uses channels 1, 2, 4 and 5. Neither 3a, nor 3b.
        """
        return get_tsm_idx(channels[:, :, 0], channels[:, :, 1],
                           channels[:, :, 4], channels[:, :, 5])


def main_klm(reader_cls, filename, start_line, end_line):
    """Generate a l1c file."""
    tic = datetime.datetime.now()
    reader = reader_cls.fromfile(filename)
    reader.save(int(start_line), int(end_line))
    LOG.info("pygac took: %s", str(datetime.datetime.now() - tic))
