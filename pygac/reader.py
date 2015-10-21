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

"""Generic reader for GAC and LAC data. Can't be used as is, has to be subclassed to add
specific read functions.
"""
import ConfigParser
import os
import logging
import datetime

import numpy as np
from pyorbital.orbital import Orbital

from pyorbital import astronomy

from pygac import CONFIG_FILE
from pygac.calibration import calibrate_solar, calibrate_thermal


# rpy values from here:http://yyy.rsmas.miami.edu/groups/rrsl/pathfinder/Processing/proc_app_a.html
rpy_coeffs = {
    'noaa7':  {'roll':  0.000,
               'pitch': 0.000,
               'yaw':   0.000,
               },
    'noaa9':  {'roll':  0.000,
               'pitch': 0.0025,
               'yaw':   0.000,
               },
    'noaa10': {'roll':  0.000,
               'pitch': 0.000,
               'yaw':   0.000,
               },
    'noaa11': {'roll':  -0.0019,
               'pitch': -0.0037,
               'yaw':   0.000,
               },
    'noaa12': {'roll':  0.000,
               'pitch': 0.000,
               'yaw':   0.000,
               },
    'noaa14': {'roll':  0.000,
               'pitch': 0.000,
               'yaw':   0.000,
               }}

LOG = logging.getLogger(__name__)


class Reader(object):
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
        self.scan_width = None
        self.scan_points = None

    def get_counts(self):
        packed_data = self.scans["sensor_data"]
        counts = np.zeros((len(self.scans), self.scan_width * 5))
        counts_nb = (self.scan_width * 5) / 3
        remainder = (self.scan_width * 5) % 3
        if remainder == 0:
            nb1 = nb2 = nb3 = counts_nb
        elif remainder == 1:
            nb1 = counts_nb + 1
            nb2 = nb3 = counts_nb
        elif remainder == 2:
            nb1 = nb2 = counts_nb + 1
            nb3 = counts_nb

        counts[:, 0::3] = ((packed_data >> 20) & 1023)[:, :nb1]
        counts[:, 1::3] = ((packed_data >> 10) & 1023)[:, :nb2]
        counts[:, 2::3] = (packed_data & 1023)[:, :nb3]
        counts = counts.reshape((-1, self.scan_width, 5))
        try:
            switch = self.get_ch3_switch()
        except AttributeError:
            return counts
        else:
            channels = np.zeros((len(self.scans), self.scan_width, 6),
                                dtype=counts.dtype)
            channels[:, :, :2] = counts[:, :, :2]
            channels[:, :, -2:] = counts[:, :, -2:]
            channels[:, :, 2][switch == 1] = counts[:, :, 2][switch == 1]
            channels[:, :, 3][switch == 0] = counts[:, :, 2][switch == 0]
            return channels

    def get_times(self):
        raise NotImplementedError

    def compute_lonlat(self, utcs=None, clock_drift_adjust=True):
        if utcs is None:
            utcs = self.get_times()

        tle1, tle2 = self.get_tle_lines()

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
            LOG.debug("Using attitude correction from file")
            rpy = np.deg2rad([self.head["constant_roll_attitude_error"] / 1e3,
                              self.head["constant_pitch_attitude_error"] / 1e3,
                              self.head["constant_yaw_attitude_error"] / 1e3])
        else:
            try:
                rpy_spacecraft = rpy_coeffs[self.spacecraft_name]
                rpy = [rpy_spacecraft['roll'], rpy_spacecraft['pitch'], rpy_spacecraft['yaw']]
                LOG.debug("Using static attitude correction")
            except KeyError:
                LOG.debug("Not applying attitude correction")
                rpy = [0, 0, 0]


        LOG.info("Using rpy: %s", str(rpy))

        from pyorbital.geoloc_instrument_definitions import avhrr_gac
        from pyorbital.geoloc import compute_pixels, get_lonlatalt
        sgeom = avhrr_gac(utcs.astype(datetime.datetime), self.scan_points, 55.385)
        s_times = sgeom.times(t)

        pixels_pos = compute_pixels((tle1, tle2), sgeom, s_times, rpy)
        pos_time = get_lonlatalt(pixels_pos, s_times)

        toc = datetime.datetime.now()

        LOG.warning("Computation of geolocation: %s", str(toc - tic))

        lons, lats = pos_time[:2]

        return lons.reshape(-1, len(self.scan_points)), lats.reshape(-1, len(self.scan_points))

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

    def adjust_clock_drift(self):
        pass

class GACReader(Reader):
    def __init__(self):
        super(GACReader, self).__init__()
        self.scan_width = 409
        self.scan_points = np.arange(3.5, 2048, 5)



class LACReader(Reader):
    def __init__(self):
        super(LACReader, self).__init__()
        self.scan_width = 2048
        self.scan_points = np.arange(24, 2048, 40)

