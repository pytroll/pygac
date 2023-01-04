#!/usr/bin/python
# Copyright (c) 2014-2019 Pygac developers
#

# Author(s):

#   Abhay Devasthale <abhay.devasthale@smhi.se>
#   Adam Dybbroe <adam.dybbroe@smhi.se>
#   Sajid Pareeth <sajid.pareeth@fmach.it>
#   Martin Raspaud <martin.raspaud@smhi.se>
#   Carlos Horn <carlos.horn@external.eumetsat.int>

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

Reads L1b GAC/LAC data from POD series of satellites (NOAA-14 and earlier).
Format specification can be found in chapters 2 & 3 of the `POD user guide`_.

.. _POD user guide:
    https://www.ncei.noaa.gov/pub/data/satellite/publications/podguides/TIROS-N%20thru%20N-14/
"""

import datetime
import logging
try:
    from enum import IntFlag
except ImportError:
    # python version < 3.6, use a simple object without nice representation
    IntFlag = object

import numpy as np

from pyorbital.geoloc_instrument_definitions import avhrr_gac
from pyorbital.geoloc import compute_pixels, get_lonlatalt

from pygac.clock_offsets_converter import get_offsets
from pygac.correct_tsm_issue import TSM_AFFECTED_INTERVALS_POD, get_tsm_idx
from pygac.reader import Reader, ReaderError, NoTLEData
from pygac.slerp import slerp
from pygac.utils import file_opener

LOG = logging.getLogger(__name__)


class POD_QualityIndicator(IntFlag):
    """Quality Indicators.

    Source:
        POD guide Table 3.1.2.1-2. Format of quality indicators.
    """

    # POD guide Table 3.1.2.1-2. Format of quality indicators.
    FATAL_FLAG = 2**31  # Data should not be used for product generation
    TIME_ERROR = 2**30  # A time sequence error was detected while Processing
    # this frame
    DATA_GAP = 2**29  # A gap precedes this frame
    DATA_JITTER = 2**28  # Resync occurred on this frame
    CALIBRATION = 2**27  # Insufficient data for calibration
    NO_EARTH_LOCATION = 2**26  # Earth location data not available
    ASCEND_DESCEND = 2**25  # AVHRR Earth location indication of Ascending (=0)
    # or descending (=1) data
    PSEUDO_NOISE = 2**24  # Pseudo Noise (P/N) occurred (=1) on the frame,
    # data not used for calibration computations
    BIT_SYNC_STATUS = 2**23  # Drop lock during frame
    SYNC_ERROR = 2**22  # Frame Sync word error greater than zero
    FRAME_SYNC_LOCK = 2**21  # Frame Sync previously dropped lock
    FLYWHEELING = 2**20  # Flywheeling detected during this frame
    BIT_SLIPPAGE = 2**19  # Bit slippage detected during this frame
    # Solar blackbody contamination indicator
    # 0 = no correction
    # 1 = solar contamination corrected
    CH_3_CONTAMINATION = 2**18  # Channel 3 solar blackbody contamination
    CH_4_CONTAMINATION = 2**17  # Channel 4 solar blackbody contamination
    CH_5_CONTAMINATION = 2**16  # Channel 5 solar blackbody contamination
    # TIP Parity
    TIP_PARITY_1 = 2**15  # In first minor frame
    TIP_PARITY_2 = 2**14  # In second minor frame
    TIP_PARITY_3 = 2**13  # In third minor frame
    TIP_PARITY_4 = 2**12  # In fourth minor frame
    TIP_PARITY_5 = 2**11  # In fifth minor frame
    # Note: Bit 10 to 8, and 1 to 0 are spare bits. 7 to 2 define
    #       "SYNC ERRORS - Number of bit errors in frame sync" (6 bit integer?)


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

    QFlag = POD_QualityIndicator
    _quality_indicators_key = "quality_indicators"

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

    def read(self, filename, fileobj=None):
        """Read the data.

        Args:
            filename (str): Path to GAC/LAC file
            fileobj: An open file object to read from. (optional)

        Returns:
            header: numpy record array
                The header metadata
            scans: numpy record array
                The scanlines

        """
        self.filename = filename
        LOG.info('Reading %s', self.filename)
        # choose the right header depending on the date
        with file_opener(fileobj or filename) as fd_:
            self.tbm_head, self.head = self.read_header(
                filename, fileobj=fd_)
            if self.tbm_head:
                tbm_offset = tbm_header.itemsize
            else:
                tbm_offset = 0
            # read scan lines until end of file
            fd_.seek(self.offset + tbm_offset, 0)
            buffer = fd_.read()
            count = self.head["number_of_scans"]
            self._read_scanlines(buffer, count)
        year, jday, _ = self.decode_timestamps(self.head["start_time"])
        start_date = (datetime.date(year, 1, 1) +
                      datetime.timedelta(days=int(jday) - 1))
        self.correct_scan_line_numbers()
        self.spacecraft_id = self.head["noaa_spacecraft_identification_code"]
        if self.spacecraft_id == 1 and start_date < datetime.date(1982, 1, 1):
            self.spacecraft_id = 25
        self.spacecraft_name = self.spacecraft_names[self.spacecraft_id]
        LOG.info(
            "Reading %s data", self.spacecrafts_orbital[self.spacecraft_id])
        return self.head, self.scans

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
        # choose the right header depending on the date
        with file_opener(fileobj or filename) as fd_:
            # read tbm_header if present
            _tbm_head, = np.frombuffer(
                fd_.read(tbm_header.itemsize),
                dtype=tbm_header, count=1)
            try:
                data_set_name = _tbm_head['data_set_name'].decode()
            except UnicodeDecodeError:
                data_set_name = '---'
            allowed_empty = (42*b'\x00' + b'  ')
            if (cls.data_set_pattern.match(data_set_name)
                    or (_tbm_head['data_set_name'] == allowed_empty)):
                tbm_head = _tbm_head.copy()
                tbm_offset = tbm_header.itemsize
            else:
                fd_.seek(0)
                tbm_head = None
                tbm_offset = 0
            # read header
            head0, = np.frombuffer(
                fd_.read(header0.itemsize),
                dtype=header0, count=1)
            year, jday, _ = cls.decode_timestamps(head0["start_time"])
            start_date = (datetime.date(year, 1, 1) +
                          datetime.timedelta(days=int(jday) - 1))
            if start_date < datetime.date(1992, 9, 8):
                header = header1
            elif start_date <= datetime.date(1994, 11, 15):
                header = header2
            else:
                header = header3
            fd_.seek(tbm_offset, 0)
            # need to copy frombuffer to have write access on head
            head, = np.frombuffer(
                fd_.read(header.itemsize),
                dtype=header, count=1).copy()
        head = cls._correct_data_set_name(head, filename)
        cls._validate_header(head)
        return tbm_head, head

    @classmethod
    def _validate_header(cls, header):
        """Check if the header belongs to this reader."""
        # call super to enter the Method Resolution Order (MRO)
        super(PODReader, cls)._validate_header(header)
        LOG.debug("validate header")
        data_set_name = header['data_set_name'].decode()
        # split header into parts
        creation_site, transfer_mode, platform_id = (
            data_set_name.split('.')[:3])
        allowed_ids = ['TN', 'NA', 'NB', 'NC', 'ND', 'NE', 'NF', 'NG',
                       'NH', 'NI', 'NJ']
        if platform_id not in allowed_ids:
            raise ReaderError('Improper platform id "%s"!' % platform_id)

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

    def _compute_missing_lonlat(self, missed_utcs):
        """Compute lon lat values using pyorbital."""
        tic = datetime.datetime.now()

        scan_rate = datetime.timedelta(milliseconds=1/self.scan_freq).total_seconds()
        sgeom = avhrr_gac(missed_utcs.astype(datetime.datetime),
                          self.scan_points, frequency=scan_rate)
        t0 = missed_utcs[0].astype(datetime.datetime)
        s_times = sgeom.times(t0)
        tle1, tle2 = self.get_tle_lines()

        rpy = self.get_attitude_coeffs()
        pixels_pos = compute_pixels((tle1, tle2), sgeom, s_times, rpy)
        pos_time = get_lonlatalt(pixels_pos, s_times)

        missed_lons, missed_lats = pos_time[:2]

        pixels_per_line = self.lats.shape[1]
        missed_lons = missed_lons.reshape(-1, pixels_per_line)
        missed_lats = missed_lats.reshape(-1, pixels_per_line)

        toc = datetime.datetime.now()
        LOG.warning("Computation of geolocation: %s", str(toc - tic))

        return missed_lons, missed_lats

    def _adjust_clock_drift(self):
        """Adjust the geolocation to compensate for the clock error."""
        tic = datetime.datetime.now()
        self.get_times()
        try:
            error_utcs, clock_error = get_offsets(self.spacecraft_name)
        except KeyError:
            LOG.info("No clock drift info available for %s",
                     self.spacecraft_name)
            return

        error_utcs = np.array(error_utcs, dtype='datetime64[ms]')
        # interpolate to get the clock offsets at the scan line utcs
        # the clock_error is given in seconds, so offsets are in seconds, too.
        offsets = np.interp(self.utcs.astype(np.uint64),
                            error_utcs.astype(np.uint64),
                            clock_error)
        LOG.info("Adjusting for clock drift of %s to %s seconds.",
                 str(min(offsets)),
                 str(max(offsets)))

        # For the interpolation of geolocations, we need to prepend/append
        # lines. The conversion from offset to line numbers is given by the
        # scan rate = 1/scan frequency.
        scan_rate = datetime.timedelta(milliseconds=1/self.scan_freq)
        offset_lines = offsets / scan_rate.total_seconds()

        # To avoid trouble with missing lines, we will construct an
        # array that covers the whole interpolation range.
        scan_lines = self.scans["scan_line_number"]
        shifted_lines = scan_lines - offset_lines
        shifted_lines_floor = np.floor(shifted_lines).astype(int)
        # compute the line range, note that the max requires a "+1"
        # to cover the interpolation range because of the floor.
        min_line = min(scan_lines.min(), shifted_lines_floor.min())
        max_line = max(scan_lines.max(), shifted_lines_floor.max()+1)
        num_lines = max_line - min_line + 1
        missed_lines = np.setdiff1d(np.arange(min_line, max_line+1), scan_lines)
        missed_utcs = ((missed_lines - scan_lines[0])*np.timedelta64(scan_rate, "ms")
                       + self.utcs[0])
        # calculate the missing geo locations
        try:
            missed_lons, missed_lats = self._compute_missing_lonlat(missed_utcs)
        except NoTLEData as err:
            LOG.warning('Cannot perform clock drift correction: %s', str(err))
            return

        # create arrays of lons and lats for interpolation. The locations
        # correspond to not yet corrected utcs, i.e. the time difference from
        # one line to the other should be equal to the scan rate.
        pixels_per_line = self.lats.shape[1]
        complete_lons = np.full((num_lines, pixels_per_line), np.nan,
                                dtype=np.float64)
        complete_lats = np.full((num_lines, pixels_per_line), np.nan,
                                dtype=np.float64)

        complete_lons[scan_lines - min_line] = self.lons
        complete_lats[scan_lines - min_line] = self.lats
        complete_lons[missed_lines - min_line] = missed_lons
        complete_lats[missed_lines - min_line] = missed_lats

        # perform the slerp interpolation to the corrected utc times
        slerp_t = shifted_lines - shifted_lines_floor  # in [0, 1)
        slerp_res = slerp(complete_lons[shifted_lines_floor - min_line, :],
                          complete_lats[shifted_lines_floor - min_line, :],
                          complete_lons[shifted_lines_floor - min_line + 1, :],
                          complete_lats[shifted_lines_floor - min_line + 1, :],
                          slerp_t[:, np.newaxis, np.newaxis])

        # set corrected values
        self.lons = slerp_res[:, :, 0]
        self.lats = slerp_res[:, :, 1]
        self.utcs -= (offsets * 1000).astype('timedelta64[ms]')

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

    @staticmethod
    def _get_ir_channels_to_calibrate():
        ir_channels_to_calibrate = [3, 4, 5]
        return ir_channels_to_calibrate

    def postproc(self, channels):
        """No POD specific postprocessing to be done."""
        pass

    def get_tsm_pixels(self, channels):
        """Determine pixels affected by the scan motor issue.

        Uses channels 1, 2, 4 and 5. Neither 3a, nor 3b.
        """
        return get_tsm_idx(channels[:, :, 0], channels[:, :, 1],
                           channels[:, :, 3], channels[:, :, 4])

    def _get_calibrated_channels_uniform_shape(self):
        """Prepare the channels as input for gac_io.save_gac."""
        _channels = self.get_calibrated_channels()
        # prepare input
        # maybe there is a better (less memory requiring) method
        shape = _channels.shape[:-1] + (6,)
        # empty_like does not know the kwarg shape in older versions
        channels = np.empty(shape, dtype=_channels.dtype)
        channels[:, :, 0] = _channels[:, :, 0]
        channels[:, :, 1] = _channels[:, :, 1]
        channels[:, :, 2] = np.nan
        channels[:, :, 3] = _channels[:, :, 2]
        channels[:, :, 4] = _channels[:, :, 3]
        channels[:, :, 5] = _channels[:, :, 4]
        return channels


def main_pod(reader_cls, filename, start_line, end_line):
    """Generate a l1c file."""
    tic = datetime.datetime.now()
    reader = reader_cls.fromfile(filename)
    reader.save(int(start_line), int(end_line))
    LOG.info("pygac took: %s", str(datetime.datetime.now() - tic))
