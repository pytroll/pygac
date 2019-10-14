#!/usr/bin/env python

# -*- coding: utf-8 -*-

# Copyright (c) 2014, 2015 Abhay Devasthale

# Author(s):

#   Abhay Devasthale <abhay.devasthale@smhi.se>
#   Adam Dybbroe <adam.dybbroe@smhi.se>
#   Martin Raspaud <martin.raspaud@smhi.se>

# This work was done in the framework of ESA-CCI-Clouds phase I


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

Format specification can be found here:
http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/podug/html/c2/sec2-0.htm
http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/podug/html/c3/sec3-1.htm

"""

from __future__ import print_function

import numpy as np
from pygac.gac_reader import GACReader, inherit_doc
import pygac.geotiepoints as gtp
from .correct_tsm_issue import TSM_AFFECTED_INTERVALS_POD
import datetime
from pygac import gac_io

import logging
LOG = logging.getLogger(__name__)

# common header
header0 = np.dtype([("noaa_spacecraft_identification_code", ">u1"),
                    ("data_type_code", ">u1"),
                    ("start_time", ">u2", (3, )),
                    ("number_of_scans", ">u2"),
                    ("end_time", ">u2", (3, ))])

# until 8 september 1992
header1 = np.dtype([("noaa_spacecraft_identification_code", ">u1"),
                    ("data_type_code", ">u1"),
                    ("start_time", ">u2", (3, )),
                    ("number_of_scans", ">u2"),
                    ("end_time", ">u2", (3, )),
                    ("processing_block_id", "S7"),
                    ("ramp_auto_calibration", ">u1"),
                    ("number_of_data_gaps", ">u2"),
                    ("dacs_quality", ">u1", (6, )),
                    ("calibration_parameter_id", ">i2"),
                    ("dacs_status", ">u1"),
                    ("spare1", ">i1", (5, )),
                    ("data_set_name", "S44"),
                    ("spare2", ">u2", (1568, ))])

header2 = np.dtype([("noaa_spacecraft_identification_code", ">u1"),
                    ("data_type_code", ">u1"),
                    ("start_time", ">u2", (3, )),
                    ("number_of_scans", ">u2"),
                    ("end_time", ">u2", (3, )),
                    ("processing_block_id", "S7"),
                    ("ramp_auto_calibration", ">u1"),
                    ("number_of_data_gaps", ">u2"),
                    ("dacs_quality", ">u1", (6, )),
                    ("calibration_parameter_id", ">i2"),
                    ("dacs_status", ">u1"),
                    ("spare1", ">i1", (5, )),
                    ("data_set_name", "S42"),
                    ("blankfill", "S2"),
                    ("year_of_epoch", ">u2"),
                    ("julian_day_of_epoch", ">u2"),
                    ("millisecond_utc_epoch_time_of_day", ">u4"),
                    # Keplerian orbital elements
                    ("semi_major_axis", ">f8"),
                    ("eccentricity", ">f8"),
                    ("inclination", ">f8"),
                    ("argument_of_perigee", ">f8"),
                    ("right_ascension", ">f8"),
                    ("mean_anomaly", ">f8"),
                    # cartesian inertial true date of elements
                    ("x_component_of_position_vector", ">f8"),
                    ("y_component_of_position_vector", ">f8"),
                    ("z_component_of_position_vector", ">f8"),
                    ("x_dot_component_of_position_vector", ">f8"),
                    ("y_dot_component_of_position_vector", ">f8"),
                    ("z_dot_component_of_position_vector", ">f8"),
                    ("spare2", ">u2", (1516, ))])

header3 = np.dtype([("noaa_spacecraft_identification_code", ">u1"),
                    ("data_type_code", ">u1"),
                    ("start_time", ">u2", (3, )),
                    ("number_of_scans", ">u2"),
                    ("end_time", ">u2", (3, )),
                    ("processing_block_id", "S7"),
                    ("ramp_auto_calibration", ">u1"),
                    ("number_of_data_gaps", ">u2"),
                    ("dacs_quality", ">u1", (6, )),
                    ("calibration_parameter_id", ">i2"),
                    ("dacs_status", ">u1"),
                    ("reserved_for_mounting_and_fixed_attitude_correction_indicator",
                     ">i1"),
                    ("nadir_earth_location_tolerance", ">i1"),
                    ("spare1", ">i1"),
                    ("start_of_data_set_year", ">u2"),
                    ("data_set_name", "S44"),
                    ("year_of_epoch", ">u2"),
                    ("julian_day_of_epoch", ">u2"),
                    ("millisecond_utc_epoch_time_of_day", ">u4"),
                    # Keplerian orbital elements
                    ("semi_major_axis", ">i4"),
                    ("eccentricity", ">i4"),
                    ("inclination", ">i4"),
                    ("argument_of_perigee", ">i4"),
                    ("right_ascension", ">i4"),
                    ("mean_anomaly", ">i4"),
                    # cartesian inertial true date of elements
                    ("x_component_of_position_vector", ">i4"),
                    ("y_component_of_position_vector", ">i4"),
                    ("z_component_of_position_vector", ">i4"),
                    ("x_dot_component_of_position_vector", ">i4"),
                    ("y_dot_component_of_position_vector", ">i4"),
                    ("z_dot_component_of_position_vector", ">i4"),
                    # future use
                    ("yaw_fixed_error_correction", ">i2"),
                    ("roll_fixed_error_correction", ">i2"),
                    ("pitch_fixed_error_correction", ">i2"),
                    ("spare2", ">u2", (1537, ))])


scanline = np.dtype([("scan_line_number", ">i2"),
                     ("time_code", ">u2", (3, )),
                     # ("time_code", ">u1", (6, )),
                     ("quality_indicators", ">u4"),
                     ("calibration_coefficients", ">i4", (10, )),
                     ("number_of_meaningful_zenith_angles_and_earth_location_appended",
                      ">u1"),
                     ("solar_zenith_angles", "i1", (51, )),
                     ("earth_location", ">i2", (102, )),
                     ("telemetry", ">u4", (35, )),
                     ("sensor_data", ">u4", (682, )),
                     ("add_on_zenith", ">u2", (10, )),
                     ("clock_drift_delta", ">u2"),
                     ("spare3", "u2", (11, ))])


@inherit_doc
class GACPODReader(GACReader):

    spacecrafts_orbital = {25: 'tiros n',
                           2: 'noaa 6',
                           4: 'noaa 7',
                           6: 'noaa 8',
                           7: 'noaa 9',
                           8: 'noaa 10',
                           1: 'noaa 11',
                           5: 'noaa 12',
                           3: 'noaa 14',
                           }
    spacecraft_names = {25: 'tirosn',
                        2: 'noaa6',
                        4: 'noaa7',
                        6: 'noaa8',
                        7: 'noaa9',
                        8: 'noaa10',
                        1: 'noaa11',
                        5: 'noaa12',
                        3: 'noaa14',
                        }

    tsm_affected_intervals = TSM_AFFECTED_INTERVALS_POD

    def correct_scan_line_numbers(self, plot=False):
        # Perform common corrections first.
        super(GACPODReader, self).correct_scan_line_numbers(plot=plot)

        # cleaning up the data
        min_scanline_number = np.amin(np.absolute(self.scans["scan_line_number"][:]))
        if self.scans["scan_line_number"][0] == self.scans["scan_line_number"][-1] + 1:
            while self.scans["scan_line_number"][0] != min_scanline_number:
                self.scans = np.roll(self.scans, -1)
        else:
            while self.scans["scan_line_number"][0] != min_scanline_number:
                self.scans = self.scans[1:]

        self.scans = self.scans[self.scans["scan_line_number"] != 0]

    def read(self, filename):
        super(GACPODReader, self).read(filename=filename)

        # choose the right header depending on the date
        with open(filename) as fd_:
            head = np.fromfile(fd_, dtype=header0, count=1)[0]
            year, jday, _ = self.decode_timestamps(head["start_time"])

            start_date = (datetime.date(year, 1, 1) + datetime.timedelta(days=int(jday) - 1))

            if start_date < datetime.date(1992, 9, 8):
                header = header1
            elif start_date <= datetime.date(1994, 11, 15):
                header = header2
            else:
                header = header3

        with open(filename) as fd_:
            self.head = np.fromfile(fd_, dtype=header, count=1)[0]
            self.scans = np.fromfile(fd_,
                                     dtype=scanline,
                                     count=self.head["number_of_scans"])
        self.correct_scan_line_numbers()
        self.spacecraft_id = self.head["noaa_spacecraft_identification_code"]
        if self.spacecraft_id == 1 and start_date < datetime.date(1982, 1, 1):
            self.spacecraft_id = 25
        self.spacecraft_name = self.spacecraft_names[self.spacecraft_id]
        LOG.info(
            "Reading %s data", self.spacecrafts_orbital[self.spacecraft_id])

        return self.head, self.scans

    def get_header_timestamp(self):
        year, jday, msec = self.decode_timestamps(self.head["start_time"])
        try:
            return self.to_datetime(self.to_datetime64(year=year, jday=jday,
                                                       msec=msec))
        except ValueError as err:
            raise ValueError('Corrupt header timestamp: {0}'.format(err))

    @staticmethod
    def decode_timestamps(encoded):
        """
        Decode timestamps.
        @return: year, day of year, milliseconds since 00:00
        """
        ndims = len(encoded.shape)
        if ndims == 1:
            # Single header timestamp
            enc0 = encoded[0]
            enc1 = encoded[1]
            enc2 = encoded[2]
        elif ndims == 2:
            # Scanline timestamps
            enc0 = encoded[:, 0]
            enc1 = encoded[:, 1]
            enc2 = encoded[:, 2]
        else:
            raise ValueError('Invalid timestamp dimension')

        year = enc0 >> 9
        year = np.where(year > 75, year + 1900, year + 2000)
        jday = (enc0 & 0x1FF)
        msec = ((np.uint32(enc1 & 2047) << 16) | (np.uint32(enc2)))

        return year, jday, msec

    def _get_times(self):
        return self.decode_timestamps(self.scans["time_code"])

    def adjust_clock_drift(self):
        """Adjust the geolocation to compensate for the clock error.

        TODO: bad things might happen when scanlines are skipped.
        """
        tic = datetime.datetime.now()
        self.get_times()
        from pygac.clock_offsets_converter import get_offsets
        try:
            offset_times, clock_error = get_offsets(self.spacecraft_name)
        except KeyError:
            LOG.info("No clock drift info available for %s",
                     self.spacecraft_name)
        else:
            offset_times = np.array(offset_times, dtype='datetime64[ms]')
            offsets = np.interp(self.utcs.astype(np.uint64),
                                offset_times.astype(np.uint64),
                                clock_error)
            LOG.info("Adjusting for clock drift of %s to %s",
                     str(min(offsets)),
                     str(max(offsets)))
            self.times = (self.utcs +
                          offsets.astype('timedelta64[s]')).astype(datetime.datetime)
            offsets *= -2

            int_offsets = np.floor(offsets).astype(np.int)

            # filling out missing geolocations with computation from pyorbital.
            line_indices = (self.scans["scan_line_number"]
                            + int_offsets)

            missed = sorted((set(line_indices) |
                             set(line_indices + 1))
                            - set(self.scans["scan_line_number"]))

            min_idx = min(line_indices)
            max_idx = max(max(line_indices),
                          max(self.scans["scan_line_number"] - min_idx)) + 1
            idx_len = max_idx - min_idx + 2

            complete_lons = np.zeros((idx_len, 409), dtype=np.float) * np.nan
            complete_lats = np.zeros((idx_len, 409), dtype=np.float) * np.nan

            complete_lons[self.scans["scan_line_number"] - min_idx] = self.lons
            complete_lats[self.scans["scan_line_number"] - min_idx] = self.lats

            missed_utcs = ((np.array(missed) - 1) * np.timedelta64(500, "ms")
                           + self.utcs[0])

            mlons, mlats = self.compute_lonlat(missed_utcs, True)

            complete_lons[missed - min_idx] = mlons
            complete_lats[missed - min_idx] = mlats

            from pygac.slerp import slerp
            off = offsets - np.floor(offsets)
            res = slerp(complete_lons[line_indices - min_idx, :],
                        complete_lats[line_indices - min_idx, :],
                        complete_lons[line_indices - min_idx + 1, :],
                        complete_lats[line_indices - min_idx + 1, :],
                        off[:, np.newaxis, np.newaxis])

            self.lons = res[:, :, 0]
            self.lats = res[:, :, 1]
            self.utcs += offsets.astype('timedelta64[s]')

        toc = datetime.datetime.now()
        LOG.debug("clock drift adjustment took %s", str(toc - tic))

    def get_lonlat(self):
        # interpolating lat-on points using PYTROLL geotiepoints
        arr_lat = self.scans["earth_location"][:, 0::2] / 128.0
        arr_lon = self.scans["earth_location"][:, 1::2] / 128.0

        self.lons, self.lats = gtp.Gac_Lat_Lon_Interpolator(arr_lon, arr_lat)
        return self.lons, self.lats

    def get_telemetry(self):
        number_of_scans = self.scans["telemetry"].shape[0]
        decode_tele = np.zeros((int(number_of_scans), 105))
        decode_tele[:, ::3] = (self.scans["telemetry"] >> 20) & 1023
        decode_tele[:, 1::3] = (self.scans["telemetry"] >> 10) & 1023
        decode_tele[:, 2::3] = self.scans["telemetry"] & 1023

        prt_counts = np.mean(decode_tele[:, 17:20], axis=1)

        # getting ICT counts

        ict_counts = np.zeros((int(number_of_scans), 3))
        ict_counts[:, 0] = np.mean(decode_tele[:, 22:50:3], axis=1)
        ict_counts[:, 1] = np.mean(decode_tele[:, 23:51:3], axis=1)
        ict_counts[:, 2] = np.mean(decode_tele[:, 24:52:3], axis=1)

        # getting space counts

        space_counts = np.zeros((int(number_of_scans), 3))
        space_counts[:, 0] = np.mean(decode_tele[:, 54:100:5], axis=1)
        space_counts[:, 1] = np.mean(decode_tele[:, 55:101:5], axis=1)
        space_counts[:, 2] = np.mean(decode_tele[:, 56:102:5], axis=1)

        return prt_counts, ict_counts, space_counts

    def get_corrupt_mask(self):

        # corrupt scanlines

        mask = ((self.scans["quality_indicators"] >> 31) |
                ((self.scans["quality_indicators"] << 4) >> 31) |
                ((self.scans["quality_indicators"] << 5) >> 31))

        number_of_scans = self.scans["telemetry"].shape[0]
        qual_flags = np.zeros((int(number_of_scans), 7))
        qual_flags[:, 0] = self.scans["scan_line_number"]
        qual_flags[:, 1] = (self.scans["quality_indicators"] >> 31)
        qual_flags[:, 2] = ((self.scans["quality_indicators"] << 4) >> 31)
        qual_flags[:, 3] = ((self.scans["quality_indicators"] << 5) >> 31)
        qual_flags[:, 4] = ((self.scans["quality_indicators"] << 13) >> 31)
        qual_flags[:, 5] = ((self.scans["quality_indicators"] << 14) >> 31)
        qual_flags[:, 6] = ((self.scans["quality_indicators"] << 15) >> 31)

        return mask.astype(bool), qual_flags


def main(filename, start_line, end_line):
    tic = datetime.datetime.now()
    reader = GACPODReader()
    reader.read(filename)
    reader.get_lonlat()
    reader.adjust_clock_drift()
    channels = reader.get_calibrated_channels()
    sat_azi, sat_zen, sun_azi, sun_zen, rel_azi = reader.get_angles()

    mask, qual_flags = reader.get_corrupt_mask()
    if (np.all(mask)):
        print("ERROR: All data is masked out. Stop processing")
        raise ValueError("All data is masked out.")

    gac_io.save_gac(reader.spacecraft_name,
                    reader.utcs,
                    reader.lats, reader.lons,
                    channels[:, :, 0], channels[:, :, 1],
                    np.ones_like(channels[:, :, 0]) * -1,
                    channels[:, :, 2],
                    channels[:, :, 3],
                    channels[:, :, 4],
                    sun_zen, sat_zen, sun_azi, sat_azi, rel_azi,
                    mask, qual_flags, start_line, end_line,
                    reader.is_tsm_affected(),
                    reader.filename,
                    reader.get_midnight_scanline(),
                    reader.get_miss_lines())
    LOG.info("pygac took: %s", str(datetime.datetime.now() - tic))


if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2], sys.argv[3])
