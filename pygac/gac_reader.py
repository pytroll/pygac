#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014, 2015 Martin Raspaud

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>

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

"""Generic reader for GAC data. Can't be used as is, has to be subclassed to add
specific read functions.
"""
import numpy as np
from pygac import CONFIG_FILE
try:
    import ConfigParser
except ImportError:
    import configparser as ConfigParser

import os
import logging
from pyorbital.orbital import Orbital
from pyorbital import astronomy
import datetime
from pygac.calibration import calibrate_solar, calibrate_thermal
from abc import ABCMeta, abstractmethod, abstractproperty
import six
import types

LOG = logging.getLogger(__name__)


class GACReader(six.with_metaclass(ABCMeta)):

    scan_freq = 2.0/1000.0
    """Scanning frequency (scanlines per millisecond)"""

    def __init__(self, tle_thresh=7):
        self.tle_thresh = tle_thresh
        self.head = None
        self.scans = None
        self.spacecraft_name = None
        self.spacecraft_id = None
        self.utcs = None
        self.lats = None
        self.lons = None
        self.times = None
        self.tle_lines = None
        self.filename = None

    @abstractmethod
    def read(self, filename):
        """Read GAC data.

        Args:
            filename (str): Specifies the GAC file to be read.
        """
        self.filename = os.path.basename(filename)
        LOG.info('Reading %s', self.filename)

    @abstractmethod
    def get_header_timestamp(self):
        """Read start timestamp from the header.

        Returns:
            datetime.datetime: Start timestamp
        """
        raise NotImplementedError

    def get_counts(self):
        packed_data = self.scans["sensor_data"]
        gac_counts = np.zeros((len(self.scans), 409 * 5))
        gac_counts[:, 0::3] = (packed_data >> 20) & 1023
        gac_counts[:, 1::3] = (packed_data >> 10) & 1023
        gac_counts[:, 2::3] = (packed_data & 1023)[:, :-1]
        gac_counts = gac_counts.reshape((-1, 409, 5))
        try:
            switch = self.get_ch3_switch()
        except AttributeError:
            return gac_counts
        else:
            channels = np.zeros((len(self.scans), 409, 6),
                                dtype=gac_counts.dtype)
            channels[:, :, :2] = gac_counts[:, :, :2]
            channels[:, :, -2:] = gac_counts[:, :, -2:]
            channels[:, :, 2][switch == 1] = gac_counts[:, :, 2][switch == 1]
            channels[:, :, 3][switch == 0] = gac_counts[:, :, 2][switch == 0]
            return channels

    @abstractmethod
    def _get_times(self):
        """Specifies how to read scanline timestamps from GAC data.

        Returns:
            int: year
            int: day of year
            int: milliseconds since 00:00
        """
        raise NotImplementedError

    def get_times(self):
        """Read scanline timestamps and try to correct invalid values.

        Note:
            Also sets self.utcs and self.times!

        Returns:
            UTC timestamps
        """
        if self.utcs is None:
            # Read timestamps
            year, jday, msec = self._get_times()

            # Correct invalid values
            year, jday, msec = self.correct_times_median(year=year, jday=jday,
                                                         msec=msec)
            self.utcs = self.to_datetime64(year=year, jday=jday, msec=msec)
            self.correct_times_thresh()

            # Convert timestamps to datetime objects
            self.times = self.to_datetime(self.utcs)

        return self.utcs

    @staticmethod
    def to_datetime64(year, jday, msec):
        """Convert timestamps to numpy.datetime64

        Args:
            year: Year
            jday: Day of the year (1-based)
            msec: Milliseconds since 00:00

        Returns:
            numpy.datetime64: Converted timestamps
        """
        return (((year - 1970).astype('datetime64[Y]')
                + (jday - 1).astype('timedelta64[D]')).astype('datetime64[ms]')
                + msec.astype('timedelta64[ms]'))

    @staticmethod
    def to_datetime(datetime64):
        """Convert numpy.datetime64 to datetime.datetime

        Args:
            datetime64 (numpy.datetime64): Numpy timestamp to be converted.

        Returns:
            datetime.datetime: Converted timestamp
        """
        return datetime64.astype(datetime.datetime)

    def lineno2msec(self, scan_line_number):
        """Compute ideal scanline timestamp based on the scanline number.

        Assumes a constant scanning frequency.

        Args:
            scan_line_number: Specifies the scanline number (1-based)

        Returns:
            Corresponding timestamps in milliseconds since 1970-01-01 00:00,
            i.e. the first scanline has timestamp 0.
        """
        return (scan_line_number - 1) / self.scan_freq

    def compute_lonlat(self, utcs=None, clock_drift_adjust=True):
        tle1, tle2 = self.get_tle_lines()

        scan_points = np.arange(3.5, 2048, 5)

        if utcs is None:
            utcs = self.get_times()

        # adjusting clock for drift
        tic = datetime.datetime.now()
        if clock_drift_adjust:
            from pygac.clock_offsets_converter import get_offsets
            try:
                offset_times, clock_error = get_offsets(self.spacecraft_name)
            except KeyError:
                LOG.info("No clock drift info available for %s",
                         self.spacecraft_name)
            else:
                offset_times = np.array(offset_times, dtype='datetime64[ms]')
                offsets = np.interp(utcs.astype(np.uint64),
                                    offset_times.astype(np.uint64),
                                    clock_error)
                utcs -= (offsets * 1000).astype('timedelta64[ms]')

        t = utcs[0].astype(datetime.datetime)

        if "constant_yaw_attitude_error" in self.head.dtype.fields:
            rpy = np.deg2rad([self.head["constant_roll_attitude_error"] / 1e3,
                              self.head["constant_pitch_attitude_error"] / 1e3,
                              self.head["constant_yaw_attitude_error"] / 1e3])
        else:
            rpy = [0, 0, 0]

        LOG.info("Using rpy: %s", str(rpy))

        from pyorbital.geoloc_instrument_definitions import avhrr_gac
        from pyorbital.geoloc import compute_pixels, get_lonlatalt
        sgeom = avhrr_gac(utcs.astype(datetime.datetime), scan_points, 55.385)
        s_times = sgeom.times(t)

        pixels_pos = compute_pixels((tle1, tle2), sgeom, s_times, rpy)
        pos_time = get_lonlatalt(pixels_pos, s_times)

        toc = datetime.datetime.now()

        LOG.warning("Computation of geolocation: %s", str(toc - tic))

        lons, lats = pos_time[:2]

        return lons.reshape(-1, 409), lats.reshape(-1, 409)

    def get_calibrated_channels(self):
        channels = self.get_counts()
        self.get_times()
        year = self.times[0].year
        delta = self.times[0].date() - datetime.date(year, 1, 1)
        jday = delta.days + 1

        # Earth-Sun distance correction factor
        corr = 1.0 - 0.0334 * np.cos(2.0 * np.pi * (jday - 2) / 365.25)

        # how many reflective channels are there ?
        tot_ref = channels.shape[2] - 3

        channels[:, :, 0:tot_ref] = calibrate_solar(channels[:, :, 0:tot_ref],
                                                    np.arange(tot_ref),
                                                    year, jday,
                                                    self.spacecraft_name,
                                                    corr)
        prt, ict, space = self.get_telemetry()
        for chan in [3, 4, 5]:
            channels[:, :, chan - 6] = calibrate_thermal(
                channels[:, :, chan - 6],
                prt,
                ict[:, chan - 3],
                space[:, chan - 3],
                self.scans["scan_line_number"],
                chan,
                self.spacecraft_name)
        return channels

    @staticmethod
    def tle2datetime64(times):
        """Convert TLE timestamps to numpy.datetime64

        Args:
           times (float): TLE timestamps as %y%j.1234, e.g. 18001.25
        """
        # Convert %y%j.12345 to %Y%j.12345 (valid for 1950-2049)
        times = np.where(times > 50000, times + 1900000, times + 2000000)

        # Convert float to datetime64
        doys = (times % 1000).astype('int') - 1
        years = (times // 1000).astype('int')
        msecs = np.rint(24 * 3600 * 1000 * (times % 1))
        times64 = (years - 1970).astype('datetime64[Y]').astype('datetime64[ms]')
        times64 += doys.astype('timedelta64[D]')
        times64 += msecs.astype('timedelta64[ms]')

        return times64

    def get_tle_file(self):
        conf = ConfigParser.ConfigParser()
        try:
            conf.read(CONFIG_FILE)
        except ConfigParser.NoSectionError:
            LOG.exception('Failed reading configuration file: %s',
                          str(CONFIG_FILE))
            raise

        values = {"satname": self.spacecraft_name, }
        options = {}
        for option, value in conf.items('tle', raw=True):
            options[option] = value
        tle_filename = os.path.join(options['tledir'],
                                    options["tlename"] % values)
        LOG.info('TLE filename = ' + str(tle_filename))
        with open(tle_filename, 'r') as fp_:
            return fp_.readlines()

    def get_tle_lines(self):
        """Find closest two line elements (TLEs) for the current orbit

        Raises:
            IndexError, if the closest TLE is more than :meth:`pygac.GACReader.tle_thresh` days apart
        """
        if self.tle_lines is not None:
            return self.tle_lines

        self.get_times()
        tle_data = self.get_tle_file()
        sdate = np.datetime64(self.times[0], '[ms]')
        dates = self.tle2datetime64(
            np.array([float(line[18:32]) for line in tle_data[::2]]))

        # Find index "iindex" such that dates[iindex-1] < sdate <= dates[iindex]
        # Notes:
        #     1. If sdate < dates[0] then iindex = 0
        #     2. If sdate > dates[-1] then iindex = len(dates), beyond the right boundary!
        iindex = np.searchsorted(dates, sdate)

        if iindex in (0, len(dates)):
            if iindex == len(dates):
                # Reset index if beyond the right boundary (see note 2. above)
                iindex -= 1
        elif abs(sdate - dates[iindex - 1]) < abs(sdate - dates[iindex]):
            # Choose the closest of the two surrounding dates
            iindex -= 1

        # Make sure the TLE we found is within the threshold
        delta_days = abs(sdate - dates[iindex]) / np.timedelta64(1, 'D')
        if delta_days > self.tle_thresh:
            raise IndexError(
                "Can't find tle data for %s within +/- %d days around %s" %
                (self.spacecraft_name, self.tle_thresh, sdate))

        if delta_days > 3:
            LOG.warning("Found TLE data for %s that is %f days appart",
                        sdate, delta_days)
        else:
            LOG.debug("Found TLE data for %s that is %f days appart",
                      sdate, delta_days)

        # Select TLE data
        tle1 = tle_data[iindex * 2]
        tle2 = tle_data[iindex * 2 + 1]
        self.tle_lines = tle1, tle2
        return tle1, tle2

    def get_angles(self):
        self.get_times()
        tle1, tle2 = self.get_tle_lines()
        orb = Orbital(self.spacecrafts_orbital[self.spacecraft_id],
                      line1=tle1, line2=tle2)

        sat_azi, sat_elev = orb.get_observer_look(self.times[:, np.newaxis],
                                                  self.lons, self.lats, 0)

        sat_zenith = 90 - sat_elev

        sun_zenith = astronomy.sun_zenith_angle(self.times[:, np.newaxis],
                                                self.lons, self.lats)

        alt, sun_azi = astronomy.get_alt_az(self.times[:, np.newaxis],
                                            self.lons, self.lats)
        del alt
        sun_azi = np.rad2deg(sun_azi)
        sun_azi = np.where(sun_azi < 0, sun_azi + 180, sun_azi - 180)

        rel_azi = abs(sat_azi - sun_azi)
        rel_azi = np.where(rel_azi > 180.0, 360.0 - rel_azi, rel_azi)

        return sat_azi, sat_zenith, sun_azi, sun_zenith, rel_azi

    def correct_times_median(self, year, jday, msec):
        """Replace invalid timestamps with statistical estimates (using median).

        Assumes that the majority of timestamps is ok.

        Args:
            year: Year
            jday: Day of the year
            msec: Milliseconds since 00:00

        Returns:
            Corrected year
            Corrected day of the year
            Corrected milliseconds
        """
        # Estimate ideal timestamps based on the scanline number. Still without
        # offset, e.g. the first scanline has timestamp 1970-01-01 00:00
        msec_lineno = self.lineno2msec(self.scans["scan_line_number"])

        jday = np.where(np.logical_or(jday < 1, jday > 366), np.median(jday), jday)
        if_wrong_jday = np.ediff1d(jday, to_begin=0)
        jday = np.where(if_wrong_jday < 0, max(jday), jday)

        if_wrong_msec = np.where(msec < 1)
        if_wrong_msec = if_wrong_msec[0]
        if len(if_wrong_msec) > 0:
            if if_wrong_msec[0] != 0:
                msec = msec[0] + msec_lineno
            else:
                msec0 = np.median(msec - msec_lineno)
                msec = msec0 + msec_lineno

        if_wrong_msec = np.ediff1d(msec, to_begin=0)
        msec = np.where(np.logical_and(np.logical_or(if_wrong_msec < -1000, if_wrong_msec > 1000), if_wrong_jday != 1),
                        msec[0] + msec_lineno, msec)

        # checking if year value is out of valid range
        if_wrong_year = np.where(
            np.logical_or(year < 1978, year > datetime.datetime.now().year))
        if_wrong_year = if_wrong_year[0]
        if len(if_wrong_year) > 0:
            # if the first scanline has valid time stamp
            if if_wrong_year[0] != 0:
                year = year[0]
                jday = jday[0]
                msec = msec[0] + msec_lineno
            # Otherwise use median time stamp
            else:
                year = np.median(year)
                jday = np.median(jday)
                msec0 = np.median(msec - msec_lineno)
                msec = msec0 + msec_lineno

        return year, jday, msec

    def correct_scan_line_numbers(self, plot=False):
        """Remove scanlines with corrupted scanline numbers

        This includes:
            - Scanline numbers outside the valide range
            - Scanline numbers deviating more than a certain threshold from the
            ideal case (1,2,3,...N)

        Example files having corrupt scanline numbers:
            - NSS.GHRR.NJ.D96144.S2000.E2148.B0720102.GC
            - NSS.GHRR.NJ.D96064.S0043.E0236.B0606162.WI
            - NSS.GHRR.NJ.D99286.S1818.E2001.B2466869.WI

        Args:
            plot (bool): If True, plot results.
        """
        along_track = np.arange(1, len(self.scans["scan_line_number"])+1)

        # Plot original scanline numbers
        if plot:
            import matplotlib.pyplot as plt
            _, (ax0, ax1) = plt.subplots(nrows=2)
            ax0.plot(along_track, self.scans["scan_line_number"], "b-",
                     label="original")

        # Remove scanlines whose scanline number is outside the valid range
        within_range = np.logical_and(self.scans["scan_line_number"] < 15000,
                                      self.scans["scan_line_number"] >= 0)
        self.scans = self.scans[within_range]

        # Remove scanlines deviating more than a certain threshold from the
        # ideal case (1,2,3,...N).
        ideal = np.arange(1, len(self.scans["scan_line_number"])+1)

        # ... Estimate possible offset (in case some scanlines are missing in
        # the beginning of the scan)
        offsets = self.scans["scan_line_number"] - ideal
        med_offset = np.median(offsets)

        # ... Compute difference to ideal case (1,2,3,...N) + offset
        diffs = np.abs(self.scans["scan_line_number"] - (ideal + med_offset))

        # ... Remove those scanlines whose difference is larger than a certain
        # threshold. For the threshold computation we only regard nonzero
        # differences.
        nz_diffs = diffs[diffs > 0]
        if len(nz_diffs) < 50:
            # Not enough differences for reliable statistics. Use fixed
            # threshold.
            thresh = 500
        else:
            mean_nz_diffs = np.mean(nz_diffs)
            std_nz_diffs = np.std(nz_diffs)
            med_nz_diffs = np.median(nz_diffs)
            mad_nz_diffs = np.median(np.abs(nz_diffs - med_nz_diffs))
            if mean_nz_diffs / float(med_nz_diffs) < 3:
                # Relatively small variation, keep (almost) everything
                thresh = mean_nz_diffs + 3*std_nz_diffs
            else:
                # Large variation, filter more agressively. Use median and
                # median absolute deviation (MAD) as they are less sensitive to
                # outliers. However, allow differences < 500 scanlines as they
                # occur quite often.
                thresh = max(500, med_nz_diffs + 3*mad_nz_diffs)
        self.scans = self.scans[diffs <= thresh]

        LOG.debug('Removed %s scanline(s) with corrupt scanline numbers',
                  str(len(along_track) - len(self.scans)))

        # Plot corrected scanline numbers
        if plot:
            along_track = along_track[within_range]
            along_track = along_track[diffs <= thresh]
            ax0.plot(along_track, self.scans["scan_line_number"], "r--",
                     label="corrected")
            ax0.set_ylabel("Scanline Number")
            ax0.set_xlabel("Along Track")
            ax0.legend(loc="best")

            ax1.plot(np.arange(len(nz_diffs)), nz_diffs)
            ax1.axhline(thresh, color="r", label="thresh={0:.2f}"
                        .format(thresh))
            ax1.set_xlabel("Index")
            ax1.set_ylabel("nonzero |n - n'|")
            ax1.legend()

            plt.tight_layout()
            plt.savefig(self.filename + "_scanline_number_correction.png",
                        bbox_inches='tight')

    def correct_times_thresh(self, max_diff_from_t0_head=6*60*1000,
                             min_frac_near_t0_head=0.01,
                             max_diff_from_ideal_t=10*1000, plot=False):
        """Correct corrupted timestamps using a threshold approach.

        The threshold approach is based on the scanline number and the header
        timestamp. It also works if the majority of scanlines has corrupted
        timestamps.

        The header timestamp is used as a guideline to estimate the offset
        between timestamps computed from the scanline number and the actual
        scanline timestamps in the data. If header timestamp and scanline
        timestamps do not match, no correction is applied.

        Once the offset has been estimated, one can calculate the ideal
        timestamps based on the scanline number. Timestamps deviating more than
        a certain threshold from the ideal timestamps are replaced by
        the ideal timestamps.

        Example files having corrupt timestamps:
            - NSS.GHRR.NA.D81193.S2329.E0116.B1061214.WI
            - NSS.GHRR.NL.D01035.S2342.E0135.B0192627.WI

        Args:
            max_diff_from_t0_head (int): Threshold for offset estimation: A
                scanline timestamp matches the header timestamp t0_head if it is
                within the interval

                    [t0_head - max_diff_from_t0_head,
                    t0_head + max_diff_from_t0_head]

                around the header timestamp.
            min_frac_near_t0_head (float): Specifies the minimum fraction of
                scanline timestamps matching the header timestamp required for
                applying the correction.
            max_diff_from_ideal_t (float): Threshold for timestamp correction:
                If a scanline timestamp deviates more than max_diff_from_ideal_t
                from the ideal timestamp, it is regarded as corrupt and will be
                replaced with the ideal timestamp.
            plot (bool): If True, plot results.
        """
        apply_corr = True
        fail_reason = ""
        dt64_msec = ">M8[ms]"

        # Check whether scanline number increases monotonically
        n = self.scans["scan_line_number"]
        if np.any(np.diff(n) < 0):
            LOG.error("Cannot perform timestamp correction. Scanline number "
                      "does not increase monotonically.")
            apply_corr = False
            fail_reason = "Scanline number jumps backwards"
            # We could already return here, but continue for plotting purpose

        # Convert time to milliseconds since 1970-01-01
        t = self.utcs.astype("i8")
        try:
            t0_head = np.array([self.get_header_timestamp().isoformat()],
                               dtype="datetime64[ms]").astype("i8")[0]
        except ValueError as err:
            LOG.error("Cannot perform timestamp correction: %s", err)
            return

        # Compute ideal timestamps based on the scanline number. Still
        # without offset, i.e. scanline 0 has timestamp 1970-01-01 00:00
        tn = self.lineno2msec(self.scans["scan_line_number"])

        # Try to determine the timestamp t0 of the first scanline. Since none
        # of the actual timestamps is trustworthy, use the header timestamp
        # as a guideline. However, the header timestamp may also be corrupted,
        # so we only apply corrections if there is a minimum fraction of
        # scanlines whose timestamps match the header timestamp.
        #
        # 1) Compute offsets between actual timestamps and idealized timestamps
        offsets = t - tn

        # 2) If the offsets of a certain minimum fraction of scanlines are
        #    within a certain interval around the header timestamp, estimate
        #    t0 by calculating the median offset among these timestamps. If not,
        #    we do not have reliable information and cannot proceed.
        near_t0_head = np.where(
            np.fabs(offsets - t0_head) <= max_diff_from_t0_head)[0]
        if near_t0_head.size / float(n.size) >= min_frac_near_t0_head:
            t0 = np.median(offsets[near_t0_head])
        else:
            LOG.error("Timestamp mismatch. Cannot perform correction.")
            fail_reason = "Timestamp mismatch"
            apply_corr = False
            t0 = 0

        # Add estimated offset to the ideal timestamps
        tn += t0

        # Replace timestamps deviating more than a certain threshold from the
        # ideal timestamp with the ideal timestamp.
        if apply_corr:
            corrupt_lines = np.where(np.fabs(t - tn) > max_diff_from_ideal_t)
            self.utcs[corrupt_lines] = tn[corrupt_lines].astype(dt64_msec)
            LOG.debug("Corrected %s timestamp(s)", str(len(corrupt_lines[0])))

        # Plot results
        if plot:
            import matplotlib.pyplot as plt
            along_track = np.arange(n.size)
            _, (ax0, ax1, ax2) = plt.subplots(nrows=3, sharex=True,
                                              figsize=(8, 10))

            ax0.plot(along_track, t, "b-", label="original")
            if apply_corr:
                ax0.plot(along_track, self.utcs, color="red", linestyle="--",
                         label="corrected")
            else:
                ax0.set_title(fail_reason)
            ax0.set_ylabel("Time [ms since 1970]")
            ax0.set_ylim((self.utcs.min()-np.timedelta64(30, "m")).astype("i8"),
                         (self.utcs.max()+np.timedelta64(30, "m")).astype("i8"))
            ax0.legend(loc="best")
            ax2.plot(along_track, n)
            ax2.set_ylabel("Scanline number")
            ax2.set_xlabel("Along Track")

            ax1.plot(along_track, offsets)
            ax1.fill_between(
                along_track,
                t0_head - np.ones(along_track.size)*max_diff_from_t0_head,
                t0_head + np.ones(along_track.size)*max_diff_from_t0_head,
                facecolor="g", alpha=0.33)
            ax1.axhline(y=t0_head, color="g", linestyle="--",
                        label="Header timestamp")
            ax1.set_ylim(t0_head-5*max_diff_from_t0_head,
                         t0_head+5*max_diff_from_t0_head)
            ax1.set_ylabel("Offset t-tn [ms]")
            ax1.legend(loc="best")

            plt.savefig(self.filename+"_timestamp_correction.png",
                        bbox_inches="tight", dpi=100)

    @abstractproperty
    def tsm_affected_intervals(self):
        """Specifies time intervals being affected by the scan motor problem.

        Returns:
            dict: Affected time intervals. A dictionary containing a list of
                (start, end) tuples for each affected platform. Both start and
                end must be datetime.datetime objects.
        """
        raise NotImplementedError

    def is_tsm_affected(self):
        """Determine whether this orbit is affected by the scan motor problem.

        Returns:
            bool: True if the orbit is affected, False otherwise.
        """
        self.get_times()
        ts = self.times[0]
        te = self.times[-1]
        try:
            for interval in self.tsm_affected_intervals[self.spacecraft_id]:
                if ts >= interval[0] and te <= interval[1]:
                    # Found a matching interval
                    return True

            # No matching interval, orbit is not affected
            return False
        except KeyError:
            # Platform is not affected at all
            return False

    def get_midnight_scanline(self):
        """Find the scanline where the UTC date increases by one day.

        Returns:
            int: The midnight scanline if it exists and is unique.
                 None, else.
        """
        self.get_times()
        d0 = np.datetime64(datetime.date(1970, 1, 1), 'D')
        days = (self.utcs.astype('datetime64[D]') - d0).astype(int)
        incr = np.where(np.diff(days) == 1)[0]
        if len(incr) != 1:
            if len(incr) > 1:
                LOG.warning('Unable to determine midnight scanline: '
                            'UTC date increases more than once. ')
            return None
        else:
            return incr[0]

    def get_miss_lines(self):
        """Find missing scanlines, i.e. scanlines which were dropped for some
        reason or were never recorded.

        Returns:
            Indices of missing scanlines
        """
        # Compare scanline number against the ideal case (1, 2, 3, ...) and
        # find the missing line numbers.
        ideal = set(range(1, self.scans['scan_line_number'][-1] + 1))
        missing = sorted(ideal.difference(set(self.scans['scan_line_number'])))
        return np.array(missing)


def inherit_doc(cls):
    """Make a class method inherit its docstring from the parent class.

    Copied from http://stackoverflow.com/a/8101598/5703449 .
    """
    for name, func in vars(cls).items():
        if isinstance(func, types.FunctionType) and not func.__doc__:
            for parent in cls.__bases__:
                parfunc = getattr(parent, name, None)
                if parfunc and getattr(parfunc, '__doc__', None):
                    func.__doc__ = parfunc.__doc__
                    break
    return cls
