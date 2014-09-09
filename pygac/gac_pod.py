#!/usr/bin/env python

# -*- coding: utf-8 -*-

# Copyright (c) 2014 Abhay Devasthale

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

import numpy as np
from pygac.gac_reader import GACReader
import pygac.geotiepoints as gtp
import datetime
from pygac import gac_io

import logging
LOG = logging.getLogger(__name__)

header = np.dtype([("noaa_spacecraft_identification_code", ">u1"),
                   ("data_type_code", ">u1"),
                   ("start_time", ">u1", (6, )),
                   ("number_of_scans", ">u2"),
                   ("end_time", ">u1", (6, )),
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
                   ("spare2", ">u2", (1537, ))])


scanline = np.dtype([("scan_line_number", ">u2"),
                     ("time_code", ">u2", (3, )),
                     #("time_code", ">u1", (6, )),
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


class PODReader(GACReader):

    spacecrafts_orbital = {4: 'noaa 7',
                           7: 'noaa 9',
                           8: 'noaa 10',
                           1: 'noaa 11',
                           5: 'noaa 12',
                           3: 'noaa 14',
                           }
    spacecraft_names = {4: 'noaa07',
                        7: 'noaa09',
                        8: 'noaa10',
                        1: 'noaa11',
                        5: 'noaa12',
                        3: 'noaa14',
                        }

    def read(self, filename):
        with open(filename) as fd_:
            self.head = np.fromfile(fd_, dtype=header, count=1)[0]
            scans = np.fromfile(fd_,
                                dtype=scanline,
                                count=self.head["number_of_scans"])

        # cleaning up the data
        if scans["scan_line_number"][0] == scans["scan_line_number"][-1] + 1:
            while scans["scan_line_number"][0] != 1:
                scans = np.roll(scans, -1)
        else:
            while scans["scan_line_number"][0] != 1:
                scans = scans[1:]

        self.scans = scans[scans["scan_line_number"] != 0]

        self.spacecraft_id = self.head["noaa_spacecraft_identification_code"]
        self.spacecraft_name = self.spacecraft_names[self.spacecraft_id]
        LOG.info(
            "Reading %s data", self.spacecrafts_orbital[self.spacecraft_id])

        return self.head, self.scans

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
            LOG.info("No clock drift info available for %s",
                     self.spacecraft_name)
        else:
            offset_times = np.array(offset_times, dtype='datetime64[ms]')
            offsets = np.interp(self.utcs.astype(np.uint64),
                                offset_times.astype(np.uint64),
                                clock_error)
            self.times = (self.utcs +
                          offsets.astype('timedelta64[s]')).astype(datetime.datetime)
            offsets *= -2
            main_offset = int(np.floor(np.mean(offsets)))
            int_offsets = (np.floor(offsets)).astype(np.int)
            # FIXME: this is wrong when scanlines are not contiguous.
            line_indices = (np.arange(len(self.utcs)) + int_offsets)
            # line_indices = (self.scans["scan_line_number"]
            #                + int_offsets)
            # print "we miss", sorted(set(self.scans["scan_line_number"]) -
            # set(line_indices))
            for i, line in enumerate(line_indices):
                if line >= 0:
                    first_index = i
                    break
            for i, line in enumerate(reversed(line_indices)):
                if line < len(line_indices) - 1:
                    last_index = -i
                    break
            offsets -= main_offset

            from pygac.slerp import slerp

            last_index = len(line_indices) + last_index
            indices = line_indices[first_index:last_index]
            res = slerp(self.lons[indices, :],
                        self.lats[indices, :],
                        self.lons[indices + 1, :],
                        self.lats[indices + 1, :],
                        offsets[first_index:last_index, np.newaxis, np.newaxis])

            self.lons = res[:, :, 0]
            self.lats = res[:, :, 1]
            print "offsets", first_index, last_index
            self.scans = self.scans[first_index:last_index]
            self.times = self.times[first_index:last_index]
            self.utcs = self.utcs[first_index:last_index]

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

        return mask.astype(bool)


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
    LOG.info("pygac took: %s", str(datetime.datetime.now() - tic))

if __name__ == "__main__":
    import sys
    main(sys.argv[1])
