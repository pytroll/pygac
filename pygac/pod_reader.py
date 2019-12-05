#!/usr/bin/python
# Copyright (c) 2014-2019 Pygac developers
#

# Author(s):

#   Abhay Devasthale <abhay.devasthale@smhi.se>
#   Adam Dybbroe <adam.dybbroe@smhi.se>
#   Sajid Pareeth <sajid.pareeth@fmach.it>
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


"""POD file reading.

Format specification can be found here:
http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/podug/html/c2/sec2-0.htm
http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/podug/html/c3/sec3-1.htm
"""

import datetime
import logging

import numpy as np

from pygac.correct_tsm_issue import TSM_AFFECTED_INTERVALS_POD, get_tsm_idx
from pygac.reader import Reader

LOG = logging.getLogger(__name__)

# common header
header0 = np.dtype([("noaa_spacecraft_identification_code", ">u1"),
                    ("data_type_code", ">u1"),
                    ("start_time", ">u2", (3, )),
                    ("number_of_scans", ">u2"),
                    ("end_time", ">u2", (3, ))])


# For L1B data before September 8, 1992
# http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/podug/html/k/app-k.htm
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
                    ("data_set_name", "S44")])

# For L1B data between October 21, 1992 to November 15, 1994
# http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/podug/html/l/app-l.htm
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
                    ("julian_year_of_epoch", ">u2"),
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
                    ("z_dot_component_of_position_vector", ">f8")])

# For L1B data post November 15, 1994
# http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/podug/html/c2/sec2-0.htm
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
                    ("pitch_fixed_error_correction", ">i2")])

# archive header
tbm_header = np.dtype([('fill', 'S30'),
                       ('data_set_name', 'S44'),
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
                       ('sensor_data_word_size', 'S2'),
                       ('fill2', 'S3')])


class PODReader(Reader):
    """The POD reader."""

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

    def correct_scan_line_numbers(self):
        """Correct the scan line numbers."""
        # Perform common corrections first.
        super(PODReader, self).correct_scan_line_numbers()

        # cleaning up the data
        min_scanline_number = np.amin(
            np.absolute(self.scans["scan_line_number"][:]))
        if self.scans["scan_line_number"][0] == self.scans["scan_line_number"][-1] + 1:
            while self.scans["scan_line_number"][0] != min_scanline_number:
                self.scans = np.roll(self.scans, -1)
        else:
            while self.scans["scan_line_number"][0] != min_scanline_number:
                self.scans = self.scans[1:]

        self.scans = self.scans[self.scans["scan_line_number"] != 0]

    def read(self, filename):
        """Read the data.

        Args:
            filename: string
                The filename to read from.

        Returns:
            header: numpy record array
                The header metadata
            scans: numpy record array
                The scanlines

        """
        super(PODReader, self).read(filename=filename)
        # choose the right header depending on the date
        with open(filename) as fd_:
            # read archive header
            self.tbm_head = np.fromfile(fd_, dtype=tbm_header, count=1)[0]
            if ((not self.tbm_head['data_set_name'].startswith(self.creation_site + b'.')) and
                    (self.tbm_head['data_set_name'] != b'\x00' * 42 + b'  ')):
                fd_.seek(0)
                self.tbm_head = None
                tbm_offset = 0
            else:
                tbm_offset = tbm_header.itemsize

            head = np.fromfile(fd_, dtype=header0, count=1)[0]
            year, jday, _ = self.decode_timestamps(head["start_time"])

            start_date = (datetime.date(year, 1, 1) +
                          datetime.timedelta(days=int(jday) - 1))

            if start_date < datetime.date(1992, 9, 8):
                header = header1
            elif start_date <= datetime.date(1994, 11, 15):
                header = header2
            else:
                header = header3

            fd_.seek(tbm_offset, 0)
            self.head = np.fromfile(fd_, dtype=header, count=1)[0]
            fd_.seek(self.offset + tbm_offset, 0)
            self.scans = np.fromfile(fd_,
                                     dtype=self.scanline_type,
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
        """Get the timestamp from the header.

        Returns:
            A datetime object containing the timestamp from the header.

        Raises:
            A ValueError if the timestamp is corrupt.

        """
        year, jday, msec = self.decode_timestamps(self.head["start_time"])
        try:
            return self.to_datetime(self.to_datetime64(year=year, jday=jday,
                                                       msec=msec))
        except ValueError as err:
            raise ValueError('Corrupt header timestamp: {0}'.format(err))

    @staticmethod
    def decode_timestamps(encoded):
        """Decode timestamps.

        Returns:
            year
            day of year
            milliseconds since 00:00

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

    def _adjust_clock_drift(self):
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

            complete_lons = np.full((idx_len, self.lats.shape[1]), np.nan,
                                    dtype=np.float)
            complete_lats = np.full((idx_len, self.lats.shape[1]), np.nan,
                                    dtype=np.float)

            complete_lons[self.scans["scan_line_number"] - min_idx] = self.lons
            complete_lats[self.scans["scan_line_number"] - min_idx] = self.lats

            missed_utcs = ((np.array(missed) - 1) * np.timedelta64(500, "ms")
                           + self.utcs[0])
            try:
                mlons, mlats = self.compute_lonlat(width=self.lats.shape[1],
                                                   utcs=missed_utcs,
                                                   clock_drift_adjust=True)
            except IndexError as err:
                LOG.warning(
                    'Cannot perform clock drift correction: %s', str(err))
                return

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

    def _get_lonlat(self):
        lats = self.scans["earth_location"][:, 0::2] / 128.0
        lons = self.scans["earth_location"][:, 1::2] / 128.0
        return lons, lats

    def get_telemetry(self):
        """Get the telemetry.

        Returns:
            prt_counts: np.array
            ict_counts: np.array
            space_counts: np.array

        """
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

    def _get_corrupt_mask(self):
        """Get mask for corrupt scanlines."""
        mask = ((self.scans["quality_indicators"] >> 31) |
                ((self.scans["quality_indicators"] << 4) >> 31) |
                ((self.scans["quality_indicators"] << 5) >> 31))
        return mask.astype(bool)

    def get_qual_flags(self):
        """Read quality flags."""
        number_of_scans = self.scans["telemetry"].shape[0]
        qual_flags = np.zeros((int(number_of_scans), 7))
        qual_flags[:, 0] = self.scans["scan_line_number"]
        qual_flags[:, 1] = (self.scans["quality_indicators"] >> 31)
        qual_flags[:, 2] = ((self.scans["quality_indicators"] << 4) >> 31)
        qual_flags[:, 3] = ((self.scans["quality_indicators"] << 5) >> 31)
        qual_flags[:, 4] = ((self.scans["quality_indicators"] << 13) >> 31)
        qual_flags[:, 5] = ((self.scans["quality_indicators"] << 14) >> 31)
        qual_flags[:, 6] = ((self.scans["quality_indicators"] << 15) >> 31)

        return qual_flags

    def postproc(self, channels):
        """No POD specific postprocessing to be done."""
        pass

    def get_tsm_pixels(self, channels):
        """Determine pixels affected by the scan motor issue.

        Uses channels 1, 2, 4 and 5. Neither 3a, nor 3b.
        """
        return get_tsm_idx(channels[:, :, 0], channels[:, :, 1],
                           channels[:, :, 3], channels[:, :, 4])


def main_pod(reader_cls, filename, start_line, end_line):
    """Generate a l1c file."""
    from pygac import gac_io
    tic = datetime.datetime.now()
    reader = reader_cls()
    reader.read(filename)
    reader.get_lonlat()
    channels = reader.get_calibrated_channels()
    sat_azi, sat_zen, sun_azi, sun_zen, rel_azi = reader.get_angles()

    qual_flags = reader.get_qual_flags()
    if (np.all(reader.mask)):
        print("ERROR: All data is masked out. Stop processing")
        raise ValueError("All data is masked out.")
    gac_io.save_gac(reader.spacecraft_name,
                    reader.utcs,
                    reader.lats, reader.lons,
                    channels[:, :, 0], channels[:, :, 1],
                    np.full_like(channels[:, :, 0], np.nan),
                    channels[:, :, 2],
                    channels[:, :, 3],
                    channels[:, :, 4],
                    sun_zen, sat_zen, sun_azi, sat_azi, rel_azi,
                    qual_flags, start_line, end_line,
                    reader.filename,
                    reader.meta_data)
    LOG.info("pygac took: %s", str(datetime.datetime.now() - tic))
