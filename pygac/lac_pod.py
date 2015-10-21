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

LOGGER = logging.getLogger(__name__)

####For L1B data before September 8, 1992
# http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/podug/html/k/app-k.htm
oldestheader = np.dtype([("noaa_spacecraft_identification_code", ">u1"),
                         ("data_type_code", ">u1"),
                         ("start_time", ">u2", (3,)),
                         ("number_of_scans", ">u2"),
                         ("end_time", ">u2", (3,)),
                         ("processing_block_id", "S7"),
                         ("ramp_auto_calibration", ">u1"),
                         ("number_of_data_gaps", ">u2"),
                         ("dacs_quality", ">u1", (6,)),
                         ("calibration_parameter_id", ">i2"),
                         ("dacs_status", ">u1"),
                         ("zero_filled", ">i1", (5,)),
                         ("data_set_name", "S44"),
                         # future use
                         ("spare1", ">u2", (3657,))])

oldestscanline = np.dtype([("scan_line_number", ">u2"),
                           ("time_code", ">u2", (3,)),
                           ("quality_indicators", ">u4"),
                           ("calibration_coefficients", ">i4", (10,)),
                           ("number_of_meaningful_zenith_angles_and_earth_location_appended",
                            ">u1"),
                           ("solar_zenith_angles", "i1", (51,)),
                           ("earth_location", ">i2", (102,)),
                           ("telemetry", ">u4", (35,)),
                           ("sensor_data", ">u4", (3414,)),
                           ("add_on_zenith", ">u2", (10,)),
                           ("spare3", "u2", (338,))])

####For L1B data between October 21, 1992 to November 15, 1994
# http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/podug/html/l/app-l.htm
oldheader = np.dtype([("noaa_spacecraft_identification_code", ">u1"),
                      ("data_type_code", ">u1"),
                      ("start_time", ">u2", (3,)),
                      ("number_of_scans", ">u2"),
                      ("end_time", ">u2", (3,)),
                      ("processing_block_id", "S7"),
                      ("ramp_auto_calibration", ">u1"),
                      ("number_of_data_gaps", ">u2"),
                      ("dacs_quality", ">u1", (6,)),
                      ("calibration_parameter_id", ">i2"),
                      ("dacs_status", ">u1"),
                      ("zero_filled", ">i1", (5,)),
                      ("data_set_name", "S42"),
                      ("blank_filled", ">i2"),
                      ("julian_year_of_epoch", ">u2"),
                      ("julian_day_of_epoch", ">u2"),
                      ("millisecond_utc_epoch_time_of_day", ">u4"),
                      # Keplerian orbital elements
                      ("semi_major_axis", ">u4", (2,)),
                      ("eccentricity", ">u4", (2,)),
                      ("inclination", ">u4", (2,)),
                      ("argument_of_perigee", ">u4", (2,)),
                      ("right_ascension", ">u4", (2,)),
                      ("mean_anomaly", ">u4", (2,)),
                      # cartesian inertial true date of elements
                      ("x_component_of_position_vector", ">u4", (2,)),
                      ("y_component_of_position_vector", ">u4", (2,)),
                      ("z_component_of_position_vector", ">u4", (2,)),
                      ("x_dot_component_of_position_vector", ">u4", (2,)),
                      ("y_dot_component_of_position_vector", ">u4", (2,)),
                      ("z_dot_component_of_position_vector", ">u4", (2,)),
                      # future use
                      ("spare2", ">u2", (3606,))])

oldscanline = np.dtype([("scan_line_number", ">u2"),
                        ("time_code", ">u2", (3,)),
                        ("quality_indicators", ">u4"),
                        ("calibration_coefficients", ">i4", (10,)),
                        ("number_of_meaningful_zenith_angles_and_earth_location_appended",
                         ">u1"),
                        ("solar_zenith_angles", "i1", (51,)),
                        ("earth_location", ">i2", (102,)),
                        ("telemetry", ">u4", (35,)),
                        ("sensor_data", ">u4", (3414,)),
                        ("add_on_zenith", ">u2", (10,)),
                        ("spare3", "u2", (338,))])

####For L1B data post November 15, 1994
# http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/podug/html/c2/sec2-0.htm
newheader = np.dtype([("noaa_spacecraft_identification_code", ">u1"),
                      ("data_type_code", ">u1"),
                      ("start_time", ">u2", (3,)),
                      ("number_of_scans", ">u2"),
                      ("end_time", ">u2", (3,)),
                      ("processing_block_id", "S7"),
                      ("ramp_auto_calibration", ">u1"),
                      ("number_of_data_gaps", ">u2"),
                      ("dacs_quality", ">u1", (6,)),
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
                      ("semi_major_axis", ">u4"),
                      ("eccentricity", ">u4"),
                      ("inclination", ">u4"),
                      ("argument_of_perigee", ">u4"),
                      ("right_ascension", ">u4"),
                      ("mean_anomaly", ">u4"),
                      # cartesian inertial true date of elements
                      ("x_component_of_position_vector", ">u4"),
                      ("y_component_of_position_vector", ">u4"),
                      ("z_component_of_position_vector", ">u4"),
                      ("x_dot_component_of_position_vector", ">u4"),
                      ("y_dot_component_of_position_vector", ">u4"),
                      ("z_dot_component_of_position_vector", ">u4"),
                      # future use
                      ("yaw_fixed_error_correction", ">u2"),
                      ("roll_fixed_error_correction", ">u2"),
                      ("pitch_fixed_error_correction", ">u2"),
                      ("spare2", ">u2", (1537,))])

newscanline = np.dtype([("scan_line_number", ">u2"),
                        ("time_code", ">u2", (3,)),
                        ("quality_indicators", ">u4"),
                        ("calibration_coefficients", ">i4", (10,)),
                        ("number_of_meaningful_zenith_angles_and_earth_location_appended",
                         ">u1"),
                        ("solar_zenith_angles", "i1", (51,)),
                        ("earth_location", ">i2", (102,)),
                        ("telemetry", ">u4", (35,)),
                        ("sensor_data", ">u4", (3414,)),
                        ("add_on_zenith", ">u2", (10,)),
                        ("clock_drift_delta", ">u2"),
                        ("spare3", "u2", (337,))])


def get_time(time_code):
    year = time_code[0] >> 9
    if year <= 75:
        year += 100
    year += 1900

    jday = time_code[0] & 0x1FF
    msec = ((np.uint32(time_code[1] & 2047) << 16) |
            (np.uint32(time_code[2])))

    the_time = datetime.datetime(year, 1, 1)
    the_time += datetime.timedelta(days=jday - 1, milliseconds=msec)
    return the_time


class PODReader(LACReader):
    instrument_ids = {4: 7,
                      7: 9,
                      8: 10,
                      1: 11,
                      5: 12,
                      3: 14,
                      }

    spacecrafts_orbital = {4: 'noaa 7',
                           7: 'noaa 9',
                           8: 'noaa 10',
                           1: 'noaa 11',
                           5: 'noaa 12',
                           3: 'noaa 14',
                           }
    spacecraft_names = {4: 'noaa7',
                        7: 'noaa9',
                        8: 'noaa10',
                        1: 'noaa11',
                        5: 'noaa12',
                        3: 'noaa14',
                        }

    def read(self, filename):
        with open(filename) as fd_:
            # get the time from the header (all of them contain start_time)
            time = get_time(np.fromfile(fd_, dtype=oldestheader, count=1)[0]["start_time"])
            fd_.seek(0, 0)

            if time > datetime.datetime(1994, 11, 15, 23, 59, 0):
                header = newheader
                scanline = newscanline
            elif time < datetime.datetime(1992, 10, 21, 0, 0, 0):
                header = oldestheader
                scanline = oldestscanline
            elif (time > datetime.datetime(1992, 10, 21, 0, 0, 0)
                  and time < datetime.datetime(1994, 11, 15, 23, 59, 0)):
                header = oldheader
                scanline = oldscanline
            self.head = np.fromfile(fd_, dtype=header, count=1)[0]
            # The value below is very important, in this case:14800, from the LAC header table last row.
            fd_.seek(14800, 0)
            scans = np.fromfile(fd_, dtype=scanline, count=self.head["number_of_scans"])

        if scans["scan_line_number"][0] == scans["scan_line_number"][-1] + 1:
            while scans["scan_line_number"][0] != 1:
                scans = np.roll(scans, -1)
        else:
            while scans["scan_line_number"][0] == 1:
                scans = scans[1:]
        self.scans = scans[scans["scan_line_number"] != 0]

        self.spacecraft_id = self.head["noaa_spacecraft_identification_code"]
        self.instrument_id = self.instrument_ids[self.spacecraft_id]
        self.spacecraft_name = self.spacecraft_names[self.spacecraft_id]

    def get_times(self):
        if self.utcs is None:
            year = self.scans["time_code"][:, 0] >> 9
            year = np.where(year > 75, year + 1900, year + 2000)
            jday = (self.scans["time_code"][:, 0] & 0x1FF)
            msec = ((np.uint32(self.scans["time_code"][:, 1] & 2047) << 16) |
                    (np.uint32(self.scans["time_code"][:, 2])))

            self.utcs = (((year - 1970).astype('datetime64[Y]')
                          + (jday - 1).astype('timedelta64[D]')).astype('datetime64[ms]')
                         + msec.astype('timedelta64[ms]'))
            self.times = self.utcs.astype(datetime.datetime)
        return self.utcs

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
            LOGGER.info("No clock drift info available for %s",
                        self.spacecraft_name)
        else:
            offset_times = np.array(offset_times, dtype='datetime64[ms]')
            offsets = np.interp(self.utcs.astype(np.uint64),
                                offset_times.astype(np.uint64),
                                clock_error)
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
                          max(line_indices - min_idx)) + 1
            idx_len = max_idx - min_idx + 2
            # import pdb;pdb.set_trace()
            complete_lons = np.zeros((idx_len, 2048), dtype=np.float) * np.nan
            complete_lats = np.zeros((idx_len, 2048), dtype=np.float) * np.nan

            complete_lons[line_indices - min_idx] = self.lons
            complete_lats[line_indices - min_idx] = self.lats
            missed_utcs = ((np.array(missed) - 1) * np.timedelta64(500, "ms")
                           + self.utcs[0])
            # calling compute_lonlat from lac_reader
            mlons, mlats = self.compute_lonlat(missed_utcs, True)
            nlons, nlats = gtp.lac_geo_interpolator(mlons, mlats)

            complete_lons[missed - min_idx] = nlons
            complete_lats[missed - min_idx] = nlats

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
        LOGGER.debug("clock drift adjustment took %s", str(toc - tic))

    def get_lonlat(self):
        arr_lon, arr_lat = LACReader.compute_lonlat(self, utcs=None, clock_drift_adjust=True)
        self.lons, self.lats = gtp.lac_geo_interpolator(arr_lon, arr_lat)
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


def main(filename):
    tic = datetime.datetime.now()
    reader = PODReader()
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
    LOGGER.info("pygac took: %s", str(datetime.datetime.now() - tic))
