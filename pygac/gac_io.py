#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2012, 2014 Abhay Devasthale

# Author(s):

#   Abhay Devasthale <abhay.devasthale@smhi.se>
#   Adam Dybbroe <adam.dybbroe@smhi.se>
#   Sara Hornquist <sara.hornquist@smhi.se>
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


import calendar
import datetime
import logging
import os
import time

import h5py
import numpy as np

from .correct_tsm_issue import flag_pixels as flag_tsm_pixels

try:
    import ConfigParser
except ImportError:
    import configparser as ConfigParser


LOG = logging.getLogger(__name__)


try:
    CONFIG_FILE = os.environ['PYGAC_CONFIG_FILE']
except KeyError:
    LOG.exception('Environment variable PYGAC_CONFIG_FILE not set!')
    raise

if not os.path.exists(CONFIG_FILE) or not os.path.isfile(CONFIG_FILE):
    raise IOError(str(CONFIG_FILE) + " pointed to by the environment " +
                  "variable PYGAC_CONFIG_FILE is not a file or does not exist!")

conf = ConfigParser.ConfigParser()
try:
    conf.read(CONFIG_FILE)
except ConfigParser.NoSectionError:
    LOG.exception('Failed reading configuration file: ' + str(CONFIG_FILE))
    raise

options = {}
for option, value in conf.items('output', raw=True):
    options[option] = value

OUTDIR = options['output_dir']
OUTPUT_FILE_PREFIX = options['output_file_prefix']

SUNSATANGLES_DIR = os.environ.get('SM_SUNSATANGLES_DIR', OUTDIR)
AVHRR_DIR = os.environ.get('SM_AVHRR_DIR', OUTDIR)
QUAL_DIR = os.environ.get('SM_AVHRR_DIR', OUTDIR)
MISSING_DATA = -32001
MISSING_DATA_LATLON = -999999


def update_start_end_line(start_line, end_line, temp_start_line, temp_end_line):
    """Update user start/end lines after data has been sliced using temporary
    start/end lines.

    Returns:
        Updated start_line, updated end_line
    """
    new_start_line = max(0, start_line - temp_start_line)
    new_end_line = min(end_line, temp_end_line) - temp_start_line
    return new_start_line, new_end_line


def save_gac(satellite_name,
             xutcs,
             lats, lons,
             ref1, ref2, ref3,
             bt3, bt4, bt5,
             sun_zen, sat_zen, sun_azi, sat_azi, rel_azi,
             mask, qual_flags, start_line, end_line, tsmcorr,
             gac_file, midnight_scanline, miss_lines, switch=None):

    along_track = lats.shape[0]
    last_scan_line_number = qual_flags[-1, 0]

    # Determine scanline range requested by user
    start_line = int(start_line)
    end_line = int(end_line)
    if end_line == 0:
        # If the user specifies 0 as the last scanline, process all scanlines
        end_line = along_track

    bt3 = np.where(np.logical_or(bt3 < 170.0, bt3 > 350.0),
                   MISSING_DATA, bt3 - 273.15)
    bt4 = np.where(np.logical_or(bt4 < 170.0, bt4 > 350.0),
                   MISSING_DATA, bt4 - 273.15)
    bt5 = np.where(np.logical_or(bt5 < 170.0, bt5 > 350.0),
                   MISSING_DATA, bt5 - 273.15)

    lats = np.where(np.logical_or(lats < -90.00, lats > 90.00),
                    MISSING_DATA_LATLON, lats)
    lons = np.where(np.logical_or(lons < -180.00, lons > 180.00),
                    MISSING_DATA_LATLON, lons)

    sat_azi -= 180.0
    rel_azi = abs(rel_azi)
    rel_azi = 180.0 - rel_azi

    for array in [bt3, bt4, bt5]:
        array[array != MISSING_DATA] = 100 * array[array != MISSING_DATA]
        array[mask] = MISSING_DATA
    for array in [ref1, ref2, ref3,
                  sun_zen, sat_zen, sun_azi, sat_azi, rel_azi]:
        array *= 100
        array[mask] = MISSING_DATA
    for array in [lats, lons]:
        array[array != MISSING_DATA_LATLON] = 1000.0 * \
            array[array != MISSING_DATA_LATLON]
        array[mask] = MISSING_DATA_LATLON

    for ref in [ref1, ref2, ref3]:
        ref[ref < 0] = MISSING_DATA

    if switch is not None:
        ref3[switch == 0] = MISSING_DATA
        bt3[switch == 1] = MISSING_DATA
        ref3[switch == 2] = MISSING_DATA
        bt3[switch == 2] = MISSING_DATA

    for array in [ref1, ref2, ref3, bt3, bt4, bt5]:
        array[np.isnan(array)] = MISSING_DATA

    # Choose new temporary start/end lines if lat/lon info is invalid
    no_wrong_lat = np.where(lats != MISSING_DATA_LATLON)
    temp_start_line = min(no_wrong_lat[0])
    temp_end_line = max(no_wrong_lat[0])
    if temp_start_line > start_line:
        LOG.info('New start_line chosen (due to invalid lat/lon info) = ' + str(temp_start_line))
    if end_line > temp_end_line:
        LOG.info('New end_line chosen (due to invalid lat/lon info) = ' + str(temp_end_line))

    # Slice data using temporary start/end lines
    ref1 = ref1[temp_start_line:temp_end_line + 1, :].copy()
    ref2 = ref2[temp_start_line:temp_end_line + 1, :].copy()
    ref3 = ref3[temp_start_line:temp_end_line + 1, :].copy()
    bt3 = bt3[temp_start_line:temp_end_line + 1, :].copy()
    bt4 = bt4[temp_start_line:temp_end_line + 1, :].copy()
    bt5 = bt5[temp_start_line:temp_end_line + 1, :].copy()
    sun_zen = sun_zen[temp_start_line:temp_end_line + 1, :].copy()
    sun_azi = sun_azi[temp_start_line:temp_end_line + 1, :].copy()
    sat_zen = sat_zen[temp_start_line:temp_end_line + 1, :].copy()
    sat_azi = sat_azi[temp_start_line:temp_end_line + 1, :].copy()
    rel_azi = rel_azi[temp_start_line:temp_end_line + 1, :].copy()
    lats = lats[temp_start_line:temp_end_line + 1, :].copy()
    lons = lons[temp_start_line:temp_end_line + 1, :].copy()
    miss_lines = np.sort(np.array(
        qual_flags[0:temp_start_line, 0].tolist() +
        miss_lines.tolist() +
        qual_flags[temp_end_line+1:, 0].tolist()
    ))
    qual_flags = qual_flags[temp_start_line:temp_end_line+1, :].copy()
    xutcs = xutcs[temp_start_line:temp_end_line+1].copy()

    # Update user start/end lines to the new slice
    start_line, end_line = update_start_end_line(
        start_line=start_line,
        end_line=end_line,
        temp_start_line=temp_start_line,
        temp_end_line=temp_end_line)

    # Reading time from the body of the gac file
    start = xutcs[start_line].astype(datetime.datetime)
    end = xutcs[end_line].astype(datetime.datetime)
    startdate = start.strftime("%Y%m%d")
    starttime = start.strftime("%H%M%S%f")[:-5]
    enddate = end.strftime("%Y%m%d")
    endtime = end.strftime("%H%M%S%f")[:-5]
    jday = int(start.strftime("%j"))

    # Earth-Sun distance correction factor
    corr = 1.0 - 0.0334 * np.cos(2.0 * np.pi * (jday - 2) / 365.25)

    # Slice scanline range requested by user
    ref1 = ref1[start_line:end_line+1, :].copy()
    ref2 = ref2[start_line:end_line+1, :].copy()
    ref3 = ref3[start_line:end_line+1, :].copy()
    bt3 = bt3[start_line:end_line+1, :].copy()
    bt4 = bt4[start_line:end_line+1, :].copy()
    bt5 = bt5[start_line:end_line+1, :].copy()
    sun_zen = sun_zen[start_line:end_line+1, :].copy()
    sun_azi = sun_azi[start_line:end_line+1, :].copy()
    sat_zen = sat_zen[start_line:end_line+1, :].copy()
    sat_azi = sat_azi[start_line:end_line+1, :].copy()
    rel_azi = rel_azi[start_line:end_line+1, :].copy()
    lats = lats[start_line:end_line+1, :].copy()
    lons = lons[start_line:end_line+1, :].copy()
    qual_flags = qual_flags[start_line:end_line+1, :].copy()
    xutcs = xutcs[start_line:end_line+1].copy()

    if midnight_scanline is not None:
        # Update midnight scanline to the final scanline range
        midnight_scanline -= (temp_start_line + start_line)

        # Set midnight scanline to None if it has been removed due to invalid
        # lat/lon info (< 0) or lies outside the user defined scanline range
        if midnight_scanline < 0 or (midnight_scanline > end_line or
                                     midnight_scanline < start_line):
            midnight_scanline = None

    # Compute total number of scanlines
    total_number_of_scan_lines = end_line - start_line + 1

    # Correct for temporary scan motor issue.
    #
    # TODO: The thresholds in tsm.flag_pixels() were derived from the final
    #       pygac output, that's why the correction is applied here. It would
    #       certainly be more consistent to apply the correction in GACReader,
    #       but that requires a new threshold analysis.
    if tsmcorr:
        LOG.info('Correcting for temporary scan motor issue')
        tic = datetime.datetime.now()
        (ref1, ref2, bt3, bt4, bt5, ref3) = flag_tsm_pixels(channel1=ref1,
                                                            channel2=ref2,
                                                            channel3b=bt3,
                                                            channel4=bt4,
                                                            channel5=bt5,
                                                            channel3a=ref3,
                                                            fillv=MISSING_DATA)
        LOG.debug('TSM correction took: %s', str(datetime.datetime.now() - tic))

    avhrrGAC_io(satellite_name, xutcs, startdate, enddate, starttime, endtime,
                lats, lons, ref1, ref2, ref3, bt3, bt4, bt5,
                sun_zen, sat_zen, sun_azi, sat_azi, rel_azi, qual_flags,
                start_line, end_line, total_number_of_scan_lines,
                last_scan_line_number, corr, gac_file, midnight_scanline,
                miss_lines)


def avhrrGAC_io(satellite_name, xutcs, startdate, enddate, starttime, endtime,
                arrLat_full, arrLon_full, ref1, ref2, ref3, bt3, bt4, bt5,
                arrSZA, arrSTZ, arrSAA, arrSTA, arrRAA, qual_flags,
                start_line, end_line, total_number_of_scan_lines,
                last_scan_line_number, corr, gac_file, midnight_scanline,
                miss_lines):
    import os

    # Calculate start and end time in sec1970
    t_obj = time.strptime(startdate + starttime[0:6], "%Y%m%d%H%M%S")
    starttime_sec1970 = calendar.timegm(t_obj)
    t_obj = time.strptime(enddate + endtime[0:6], "%Y%m%d%H%M%S")
    endtime_sec1970 = calendar.timegm(t_obj)

    LOG.info('Output file prefix = ' + str(OUTPUT_FILE_PREFIX))
    LOG.info('AVHRR data will be written to ' + str(AVHRR_DIR))
    ofn = os.path.join(AVHRR_DIR, (OUTPUT_FILE_PREFIX + '_avhrr_' +
                                   satellite_name + '_99999_' +
                                   startdate + 'T' + starttime + 'Z_' +
                                   enddate + 'T' + endtime + 'Z.h5'))

    LOG.info('Filename: ' + str(os.path.basename(ofn)))

    fout = h5py.File(ofn, "w")

    dset1 = fout.create_dataset("/image1/data", dtype='int16', data=ref1)
    dset2 = fout.create_dataset("/image2/data", dtype='int16', data=ref2)
    dset3 = fout.create_dataset("/image3/data", dtype='int16', data=bt3)
    dset4 = fout.create_dataset("/image4/data", dtype='int16', data=bt4)
    dset5 = fout.create_dataset("/image5/data", dtype='int16', data=bt5)
    dset6 = fout.create_dataset("/image6/data", dtype='int16', data=ref3)
    dset7 = fout.create_dataset("/where/lat/data", dtype='int32',
                                data=arrLat_full)
    dset8 = fout.create_dataset("/where/lon/data", dtype='int32',
                                data=arrLon_full)
    del dset8
    channellist = []
    channellist.append("channel1".encode('utf8'))
    channellist.append("channel2".encode('utf8'))
    channellist.append("channel3b".encode('utf8'))
    channellist.append("channel4".encode('utf8'))
    channellist.append("channel5".encode('utf8'))
    channellist.append("channel3a".encode('utf8'))
    dset10 = fout.create_dataset("/how/channel_list",
                                 data=channellist)
    del dset10
    # Attributes directly on highest level groups
    g1 = fout.require_group("/image1")
    g2 = fout.require_group("/image2")
    g3 = fout.require_group("/image3")
    g4 = fout.require_group("/image4")
    g5 = fout.require_group("/image5")
    g6 = fout.require_group("/image6")
    g7 = fout.require_group("/where")

    g1.attrs["channel"] = np.string_("1")
    g1.attrs["description"] = np.string_("AVHRR ch1")
    g2.attrs["channel"] = np.string_("2")
    g2.attrs["description"] = np.string_("AVHRR ch2")
    g3.attrs["channel"] = np.string_("3b")
    g3.attrs["description"] = np.string_("AVHRR ch3b")
    g4.attrs["channel"] = np.string_("4")
    g4.attrs["description"] = np.string_("AVHRR ch4")
    g5.attrs["channel"] = np.string_("5")
    g5.attrs["description"] = np.string_("AVHRR ch5")
    g6.attrs["channel"] = np.string_("3a")
    g6.attrs["description"] = np.string_("AVHRR ch3a")
    g7.attrs["num_of_pixels"] = np.int32(arrSZA.shape[1])
    g7.attrs["num_of_lines"] = np.int32(arrSZA.shape[0])
    g7.attrs["xscale"] = np.float32(0.0)
    g7.attrs["yscale"] = np.float32(0.0)
    g7.attrs["start_line"] = start_line
    g7.attrs["end_line"] = end_line

    # Attributes in the 'what' groups
    g1 = fout.create_group("/image1/what")
    g2 = fout.create_group("/image2/what")
    g3 = fout.create_group("/image3/what")
    g4 = fout.create_group("/image4/what")
    g5 = fout.create_group("/image5/what")
    g6 = fout.create_group("/image6/what")
    g7 = fout.create_group("/where/lat/what")
    g8 = fout.create_group("/where/lon/what")
    g9 = fout.create_group("/what")

    g1.attrs["product"] = np.string_("SATCH")
    g1.attrs["quantity"] = np.string_("REFL")
    g1.attrs["dataset_name"] = np.string_('Channel 1 reflectance')
    g1.attrs["units"] = np.string_('%')
    g1.attrs["gain"] = np.float32(0.01)
    g1.attrs["offset"] = np.float32(0.0)
    g1.attrs["missingdata"] = np.int32(MISSING_DATA)
    g1.attrs["nodata"] = np.int32(MISSING_DATA)
    g1.attrs["starttime"] = np.string_(starttime[0:6])
    g1.attrs["endtime"] = np.string_(endtime[0:6])
    g1.attrs["startdate"] = np.string_(startdate)
    g1.attrs["enddate"] = np.string_(enddate)

    g2.attrs["product"] = np.string_("SATCH")
    g2.attrs["quantity"] = np.string_("REFL")
    g2.attrs["dataset_name"] = np.string_('Channel 2 reflectance')
    g2.attrs["units"] = np.string_('%')
    g2.attrs["gain"] = np.float32(0.01)
    g2.attrs["offset"] = np.float32(0.0)
    g2.attrs["missingdata"] = np.int32(MISSING_DATA)
    g2.attrs["nodata"] = np.int32(MISSING_DATA)
    g2.attrs["starttime"] = np.string_(starttime[0:6])
    g2.attrs["endtime"] = np.string_(endtime[0:6])
    g2.attrs["startdate"] = np.string_(startdate)
    g2.attrs["enddate"] = np.string_(enddate)

    g6.attrs["product"] = np.string_("SATCH")
    g6.attrs["quantity"] = np.string_("REFL")
    g6.attrs["dataset_name"] = np.string_('Channel 3a reflectance')
    g6.attrs["units"] = np.string_('%')
    g6.attrs["gain"] = np.float32(0.01)
    g6.attrs["offset"] = np.float32(0.0)
    g6.attrs["missingdata"] = np.int32(MISSING_DATA)
    g6.attrs["nodata"] = np.int32(MISSING_DATA)
    g6.attrs["starttime"] = np.string_(starttime[0:6])
    g6.attrs["endtime"] = np.string_(endtime[0:6])
    g6.attrs["startdate"] = np.string_(startdate)
    g6.attrs["enddate"] = np.string_(enddate)

    g3.attrs["product"] = np.string_("SATCH")
    g3.attrs["quantity"] = np.string_("TB")
    g3.attrs["dataset_name"] = np.string_('Channel 3b brightness temperature')
    g3.attrs["units"] = np.string_('K')
    g3.attrs["gain"] = np.float32(0.01)
    g3.attrs["offset"] = np.float32(273.15)
    g3.attrs["missingdata"] = np.int32(MISSING_DATA)
    g3.attrs["nodata"] = np.int32(MISSING_DATA)
    g3.attrs["starttime"] = np.string_(starttime[0:6])
    g3.attrs["endtime"] = np.string_(endtime[0:6])
    g3.attrs["startdate"] = np.string_(startdate)
    g3.attrs["enddate"] = np.string_(enddate)

    g4.attrs["product"] = np.string_("SATCH")
    g4.attrs["quantity"] = np.string_("TB")
    g4.attrs["dataset_name"] = np.string_('Channel 4 brightness temperature')
    g4.attrs["units"] = np.string_('K')
    g4.attrs["gain"] = np.float32(0.01)
    g4.attrs["offset"] = np.float32(273.15)
    g4.attrs["missingdata"] = np.int32(MISSING_DATA)
    g4.attrs["nodata"] = np.int32(MISSING_DATA)
    g4.attrs["starttime"] = np.string_(starttime[0:6])
    g4.attrs["endtime"] = np.string_(endtime[0:6])
    g4.attrs["startdate"] = np.string_(startdate)
    g4.attrs["enddate"] = np.string_(enddate)

    g5.attrs["product"] = np.string_("SATCH")
    g5.attrs["quantity"] = np.string_("TB")
    g5.attrs["dataset_name"] = np.string_('Channel 5 brightness temperature')
    g5.attrs["units"] = np.string_('K')
    g5.attrs["gain"] = np.float32(0.01)
    g5.attrs["offset"] = np.float32(273.15)
    g5.attrs["missingdata"] = np.int32(MISSING_DATA)
    g5.attrs["nodata"] = np.int32(MISSING_DATA)
    g5.attrs["starttime"] = np.string_(starttime[0:6])
    g5.attrs["endtime"] = np.string_(endtime[0:6])
    g5.attrs["startdate"] = np.string_(startdate)
    g5.attrs["enddate"] = np.string_(enddate)

    g7.attrs["dataset_name"] = np.string_('Latitude')
    g7.attrs["units"] = np.string_('Deg')
    g7.attrs["gain"] = np.float32(0.0010)
    g7.attrs["offset"] = np.float32(0.0)
    g7.attrs["missingdata"] = np.int32(MISSING_DATA_LATLON)
    g7.attrs["nodata"] = np.int32(MISSING_DATA_LATLON)
    g7.attrs["starttime"] = np.string_(starttime[0:6])
    g7.attrs["endtime"] = np.string_(endtime[0:6])
    g7.attrs["startdate"] = np.string_(startdate)
    g7.attrs["enddate"] = np.string_(enddate)

    g8.attrs["dataset_name"] = np.string_('Longitude')
    g8.attrs["units"] = np.string_('Deg')
    g8.attrs["gain"] = np.float32(0.0010)
    g8.attrs["offset"] = np.float32(0.0)
    g8.attrs["missingdata"] = np.int32(MISSING_DATA_LATLON)
    g8.attrs["nodata"] = np.int32(MISSING_DATA_LATLON)
    g8.attrs["starttime"] = np.string_(starttime[0:6])
    g8.attrs["endtime"] = np.string_(endtime[0:6])
    g8.attrs["startdate"] = np.string_(startdate)
    g8.attrs["enddate"] = np.string_(enddate)

    g9.attrs["object"] = np.string_("SATP")
    g9.attrs["sets"] = np.int32(len(channellist))
    g9.attrs["version"] = np.string_("H5rad ?.?")
    g9.attrs["date"] = np.string_(startdate)
    g9.attrs["time"] = np.string_(starttime[0:6])

    # Attributes in the 'how' groups
    g1 = fout.create_group("/image1/how")
    g2 = fout.create_group("/image2/how")
    g3 = fout.create_group("/image3/how")
    g4 = fout.create_group("/image4/how")
    g5 = fout.create_group("/image5/how")
    g6 = fout.create_group("/image6/how")
    g10 = fout.require_group("/how")

    # SHq: Is the sun_earth_distance correction applied?
    g1.attrs["sun_earth_distance_correction_applied"] = np.string_("TRUE")
    g1.attrs["sun_earth_distance_correction_factor"] = corr
    g2.attrs["sun_earth_distance_correction_applied"] = np.string_("TRUE")
    g2.attrs["sun_earth_distance_correction_factor"] = corr
    # No attributes on 'how' for image3,4,5
    g6.attrs["sun_earth_distance_correction_applied"] = np.string_("TRUE")
    g6.attrs["sun_earth_distance_correction_factor"] = corr

    # We do not know much about how; mostly use no-data
    g10.attrs["yaw_error"] = 0.0
    g10.attrs["roll_error"] = 0.0
    g10.attrs["pitch_error"] = 0.0
    g10.attrs["startepochs"] = starttime_sec1970
    g10.attrs["endepochs"] = endtime_sec1970
    g10.attrs["platform"] = np.string_(satellite_name)
    g10.attrs["instrument"] = np.string_("avhrr")
    g10.attrs["orbit_number"] = np.int32(99999)
    g10.attrs["gac_file"] = np.string_(gac_file)
    g10.attrs["software"] = np.string_("pyGAC")
    g10.attrs["version"] = np.string_("1.0")

    fout.close()

    LOG.info('Sun and Satellite viewing angles will be ' +
             'written to ' + str(SUNSATANGLES_DIR))
    ofn = os.path.join(SUNSATANGLES_DIR,
                       (OUTPUT_FILE_PREFIX + '_sunsatangles_' +
                        satellite_name + '_99999_' + startdate +
                        'T' + starttime + 'Z_' +
                        enddate + 'T' + endtime + 'Z.h5'))

    LOG.info('Filename: ' + str(os.path.basename(ofn)))
    fout = h5py.File(ofn, "w")

    dset1 = fout.create_dataset("/image1/data", dtype='int16', data=arrSZA)
    dset2 = fout.create_dataset("/image2/data", dtype='int16', data=arrSTZ)
    dset3 = fout.create_dataset("/image3/data", dtype='int16', data=arrRAA)
    dset4 = fout.create_dataset("/image4/data", dtype='int16', data=arrSAA)
    dset5 = fout.create_dataset("/image5/data", dtype='int16', data=arrSTA)
    dset6 = fout.create_dataset("/where/lat/data", dtype='int32',
                                data=arrLat_full)
    dset7 = fout.create_dataset("/where/lon/data", dtype='int32',
                                data=arrLon_full)

    del dset4, dset5, dset6, dset7

    # Attributes directly on highest level groups
    g1 = fout.require_group("/image1")
    g2 = fout.require_group("/image2")
    g3 = fout.require_group("/image3")
    g4 = fout.require_group("/image4")
    g5 = fout.require_group("/image5")
    g6 = fout.require_group("/where")

    g1.attrs["description"] = np.string_('Solar zenith angle')
    g2.attrs["description"] = np.string_('Satellite zenith angle')
    g3.attrs["description"] = np.string_(
        'Relative satellite-sun azimuth angle')
    g4.attrs["description"] = np.string_('Solar azimuth angle')
    g5.attrs["description"] = np.string_('Satellite azimuth angle')
    g6.attrs["num_of_pixels"] = np.int32(arrSZA.shape[1])
    g6.attrs["num_of_lines"] = np.int32(arrSZA.shape[0])
    g6.attrs["xscale"] = np.float32(0.0)
    g6.attrs["yscale"] = np.float32(0.0)
    g6.attrs["start_line"] = start_line
    g6.attrs["end_line"] = end_line

    # Attributes in the 'what' groups + 'how'
    g1 = fout.create_group("/image1/what")
    g2 = fout.create_group("/image2/what")
    g3 = fout.create_group("/image3/what")
    g4 = fout.create_group("/image4/what")
    g5 = fout.create_group("/image5/what")
    g6 = fout.create_group("/where/lat/what")
    g7 = fout.create_group("/where/lon/what")
    g8 = fout.create_group("/what")
    g9 = fout.create_group("/how")

    g1.attrs["product"] = np.string_("SUNZ")
    g1.attrs["quantity"] = np.string_("DEG")
    g1.attrs["dataset_name"] = np.string_('Solar zenith angle')
    g1.attrs["units"] = np.string_('Deg')
    g1.attrs["gain"] = np.float32(0.01)
    g1.attrs["offset"] = np.float32(0.0)
    g1.attrs["missingdata"] = np.int32(MISSING_DATA)
    g1.attrs["nodata"] = np.int32(MISSING_DATA)
    g1.attrs["starttime"] = np.string_(starttime[0:6])
    g1.attrs["endtime"] = np.string_(endtime[0:6])
    g1.attrs["startdate"] = np.string_(startdate)
    g1.attrs["enddate"] = np.string_(enddate)

    g2.attrs["product"] = np.string_("SATZ")
    g2.attrs["quantity"] = np.string_("DEG")
    g2.attrs["dataset_name"] = np.string_('Satellite zenith angle')
    g2.attrs["units"] = np.string_('Deg')
    g2.attrs["gain"] = np.float32(0.01)
    g2.attrs["offset"] = np.float32(0.0)
    g2.attrs["missingdata"] = np.int32(MISSING_DATA)
    g2.attrs["nodata"] = np.int32(MISSING_DATA)
    g2.attrs["starttime"] = np.string_(starttime[0:6])
    g2.attrs["endtime"] = np.string_(endtime[0:6])
    g2.attrs["startdate"] = np.string_(startdate)
    g2.attrs["enddate"] = np.string_(enddate)

    g3.attrs["product"] = np.string_("SSAZD")
    g3.attrs["quantity"] = np.string_("DEG")
    g3.attrs["dataset_name"] = np.string_(
        'Relative satellite-sun azimuth angle')
    g3.attrs["units"] = np.string_('Deg')
    g3.attrs["gain"] = np.float32(0.01)
    g3.attrs["offset"] = np.float32(0.0)
    g3.attrs["missingdata"] = np.int32(MISSING_DATA)
    g3.attrs["nodata"] = np.int32(MISSING_DATA)
    g3.attrs["starttime"] = np.string_(starttime[0:6])
    g3.attrs["endtime"] = np.string_(endtime[0:6])
    g3.attrs["startdate"] = np.string_(startdate)
    g3.attrs["enddate"] = np.string_(enddate)

    g4.attrs["product"] = np.string_("SUNA")
    g4.attrs["quantity"] = np.string_("DEG")
    g4.attrs["dataset_name"] = np.string_('Solar azimuth angle')
    g4.attrs["units"] = np.string_('Deg')
    g4.attrs["gain"] = np.float32(0.01)
    g4.attrs["offset"] = np.float32(180.0)
    g4.attrs["missingdata"] = np.int32(MISSING_DATA)
    g4.attrs["nodata"] = np.int32(MISSING_DATA)
    g4.attrs["starttime"] = np.string_(starttime[0:6])
    g4.attrs["endtime"] = np.string_(endtime[0:6])
    g4.attrs["startdate"] = np.string_(startdate)
    g4.attrs["enddate"] = np.string_(enddate)

    g5.attrs["product"] = np.string_("SATA")
    g5.attrs["quantity"] = np.string_("DEG")
    g5.attrs["dataset_name"] = np.string_('Satellite azimuth angle')
    g5.attrs["units"] = np.string_('Deg')
    g5.attrs["gain"] = np.float32(0.01)
    g5.attrs["offset"] = np.float32(180.0)
    g5.attrs["missingdata"] = np.int32(MISSING_DATA)
    g5.attrs["nodata"] = np.int32(MISSING_DATA)
    g5.attrs["starttime"] = np.string_(starttime[0:6])
    g5.attrs["endtime"] = np.string_(endtime[0:6])
    g5.attrs["startdate"] = np.string_(startdate)
    g5.attrs["enddate"] = np.string_(enddate)

    g6.attrs["dataset_name"] = np.string_('Latitude')
    g6.attrs["units"] = np.string_('Deg')
    g6.attrs["gain"] = np.float32(0.0010)
    g6.attrs["offset"] = np.float32(0.0)
    g6.attrs["missingdata"] = np.int32(MISSING_DATA_LATLON)
    g6.attrs["nodata"] = np.int32(MISSING_DATA_LATLON)
    g6.attrs["starttime"] = np.string_(starttime[0:6])
    g6.attrs["endtime"] = np.string_(endtime[0:6])
    g6.attrs["startdate"] = np.string_(startdate)
    g6.attrs["enddate"] = np.string_(enddate)

    g7.attrs["dataset_name"] = np.string_('Longitude')
    g7.attrs["units"] = np.string_('Deg')
    g7.attrs["gain"] = np.float32(0.0010)
    g7.attrs["offset"] = np.float32(0.0)
    g7.attrs["missingdata"] = np.int32(MISSING_DATA_LATLON)
    g7.attrs["nodata"] = np.int32(MISSING_DATA_LATLON)
    g7.attrs["starttime"] = np.string_(starttime[0:6])
    g7.attrs["endtime"] = np.string_(endtime[0:6])
    g7.attrs["startdate"] = np.string_(startdate)
    g7.attrs["enddate"] = np.string_(enddate)

    g8.attrs["object"] = np.string_("SATP")
    g8.attrs["sets"] = np.int32(5)
    g8.attrs["version"] = np.string_("H5rad ?.?")
    g8.attrs["date"] = np.string_(startdate)
    g8.attrs["time"] = np.string_(starttime[0:6])

    # We do not know much about how; mostly use no-data
    g9.attrs["yaw_error"] = 0.0
    g9.attrs["roll_error"] = 0.0
    g9.attrs["pitch_error"] = 0.0
    g9.attrs["startepochs"] = starttime_sec1970
    g9.attrs["endepochs"] = endtime_sec1970
    g9.attrs["platform"] = np.string_(satellite_name)
    g9.attrs["instrument"] = np.string_("avhrr")
    g9.attrs["orbit_number"] = np.int32(99999)
    g9.attrs["gac_file"] = np.string_(gac_file)
    g9.attrs["software"] = np.string_("pyGAC")
    g9.attrs["version"] = np.string_("1.0")

    fout.close()

    LOG.info('Quality flags will be ' +
             'written to ' + str(QUAL_DIR))
    ofn = os.path.join(QUAL_DIR,
                       (OUTPUT_FILE_PREFIX + '_qualflags_' +
                        satellite_name + '_99999_' + startdate +
                        'T' + starttime + 'Z_' +
                        enddate + 'T' + endtime + 'Z.h5'))

    LOG.info('Filename: ' + str(os.path.basename(ofn)))
    fout = h5py.File(ofn, "w")

    g1 = fout.require_group("/qual_flags")
    dset1 = g1.create_dataset("data", dtype='int16', data=qual_flags)
    del dset1

    g1.attrs["product"] = np.string_("QFLAG")
    g1.attrs["quantity"] = np.string_("INT")
    g1.attrs["dataset_name"] = np.string_('Scanline quality flags')
    g1.attrs["units"] = np.string_('None')
    g1.attrs["gain"] = np.int32(1)
    g1.attrs["offset"] = np.int32(0)
    g1.attrs["missingdata"] = np.int32(MISSING_DATA)
    g1.attrs["nodata"] = np.int32(MISSING_DATA)
    g1.attrs["starttime"] = np.string_(starttime[0:6])
    g1.attrs["endtime"] = np.string_(endtime[0:6])
    g1.attrs["startdate"] = np.string_(startdate)
    g1.attrs["enddate"] = np.string_(enddate)
    g1.attrs["gac_file"] = np.string_(gac_file)
    g1.attrs["total_number_of_data_records"] = total_number_of_scan_lines
    g1.attrs["last_scan_line_number"] = last_scan_line_number

    g2 = fout.require_group("/ancillary")
    dset2 = g2.create_dataset("missing_scanlines", dtype='int16',
                              data=miss_lines)
    del dset2
    dset3 = g2.create_dataset("scanline_timestamps", dtype='int64',
                              data=xutcs.astype('int64'))
    dset3.attrs['units'] = 'Milliseconds since 1970-01-01 00:00:00 UTC'
    dset3.attrs['calendar'] = 'standard'
    g2.attrs["midnight_scanline"] = np.string_(midnight_scanline)

    fout.close()
