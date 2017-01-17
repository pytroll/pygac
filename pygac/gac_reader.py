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
import ConfigParser
import os
import logging
from pyorbital.orbital import Orbital
from pyorbital import astronomy
import datetime
from pygac.gac_calibration import calibrate_solar, calibrate_thermal
from abc import ABCMeta, abstractmethod

LOG = logging.getLogger(__name__)


class GACReader(object):

    __metaclass__ = ABCMeta

    def __init__(self):
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
        self.filename = os.path.basename(filename)
        LOG.info('Reading ' + self.filename)

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

    def get_times(self):
        raise NotImplementedError

    @staticmethod
    def to_datetime64(year, jday, msec):
        """Convert day, day of year and milliseconds since 00:00 to
        numpy.datetime64"""
        return (((year - 1970).astype('datetime64[Y]')
                + (jday - 1).astype('timedelta64[D]')).astype('datetime64[ms]')
                + msec.astype('timedelta64[ms]'))

    @staticmethod
    def to_datetime(datetime64):
        """Convert numpy.datetime64 to datetime.datetime"""
        return datetime64.astype(datetime.datetime)

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
        if self.tle_lines is not None:
            return self.tle_lines
        tle_data = self.get_tle_file()
        tm = self.times[0]
        sdate = (int(tm.strftime("%Y%j")) +
                 (tm.hour +
                  (tm.minute +
                   (tm.second +
                    tm.microsecond / 1000000.0) / 60.0) / 60.0) / 24.0)

        dates = np.array([float(line[18:32]) for line in tle_data[::2]])
        dates = np.where(dates > 50000, dates + 1900000, dates + 2000000)

        iindex = np.searchsorted(dates, sdate)

        if ((iindex == 0 and abs(sdate - dates[0]) > 7) or
                (iindex == len(dates) - 1 and abs(sdate - dates[-1]) > 7)):
            raise IndexError(
                "Can't find tle data for %s on the %s" %
                (self.spacecraft_name, tm.isoformat()))

        if abs(sdate - dates[iindex - 1]) < abs(sdate - dates[iindex]):
            iindex -= 1

        if abs(sdate - dates[iindex]) > 3:
            LOG.warning("Found TLE data for %f that is %f days appart",
                        sdate, abs(sdate - dates[iindex]))
        else:
            LOG.debug("Found TLE data for %f that is %f days appart",
                      sdate, abs(sdate - dates[iindex]))

        tle1 = tle_data[iindex * 2]
        tle2 = tle_data[iindex * 2 + 1]
        self.tle_lines = tle1, tle2
        return tle1, tle2

    def get_angles(self):
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

    @abstractmethod
    def get_header_timestamp(self):
        """
        Read start timestamp from the header.
        @rtype: datetime.datetime
        """
        raise NotImplementedError

    def correct_scan_line_numbers(self, plot=False):
        """
        Remove scanlines with corrupt scanline numbers, i.e.

            - Scanline numbers outside the valide range
            - Scanline numbers deviating more than a certain threshold from the
            ideal case (1,2,3,...N)

        Example files having corrupt scanline numbers:
            - NSS.GHRR.NJ.D96144.S2000.E2148.B0720102.GC
            - NSS.GHRR.NJ.D96064.S0043.E0236.B0606162.WI
            - NSS.GHRR.NJ.D99286.S1818.E2001.B2466869.WI

        @param plot: If True, plot results.
        """
        along_track = np.arange(1, len(self.scans["scan_line_number"])+1)

        # Plot original scanline numbers
        if plot:
            import matplotlib.pyplot as plt
            fig, (ax0, ax1) = plt.subplots(nrows=2)
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

        LOG.debug('Removed {0} scanline(s) with corrupt scanline numbers'
                  .format(len(along_track) - len(self.scans)))

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

    def correct_utcs(self, max_diff_from_t0_head=6*60*1000,
                     min_frac_near_t0_head=0.01, max_diff_from_ideal_t=10*1000,
                     plot=False):
        """
        Correct corrupt scanline timestamps using the scanline number and an
        approximate scanning rate of 2 lines per second.

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

        @param max_diff_from_t0_head: Threshold for offset estimation: A
        scanline timestamp matches the header timestamp t0_head if it is
        within the interval

          [t0_head - max_diff_from_t0_head, t0_head + max_diff_from_t0_head]

        around the header timestamp.
        @param min_frac_near_t0_head: Specifies the minimum fraction of
        scanline timestamps matching the header timestamp required for
        applying the correction.
        @param max_diff_from_ideal_t: Threshold for timestamp correction: If
        a scanline timestamp deviates more than max_diff_from_ideal_t from
        the ideal timestamp, it is regarded as corrupt and will be replaced with
        the ideal timestamp.
        @param plot: If True, plot results.
        """
        apply_corr = True
        fail_reason = ""
        scanning_rate = 2.0 / 1000.0  # scanlines per millisecond
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
            LOG.error("Cannot perform timestamp correction: {0}".format(err))
            return

        # Compute ideal timestamps based on the scanline number. Still
        # without offset, i.e. scanline 0 has timestamp 1970-01-01 00:00
        tn = (self.scans["scan_line_number"] - 1) / scanning_rate

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
            LOG.debug("Corrected {0} timestamp(s)".format(
                      len(corrupt_lines[0])))

        # Plot results
        if plot:
            import matplotlib.pyplot as plt
            along_track = np.arange(n.size)
            fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, sharex=True,
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
