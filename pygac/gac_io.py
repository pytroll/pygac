#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2012, 2014 Abhay Devasthale

# Author(s):

#   Abhay Devasthale <abhay.devasthale@smhi.se>
#   Adam Dybbroe <adam.dybbroe@smhi.se>
#   Sara Hornquist <sara.hornquist@smhi.se>
#   Martin Raspaud <martin.raspaud@smhi.se>
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


import calendar
import datetime
import logging
import os
import time

import h5py
import numpy as np

from pygac.configuration import get_config
from pygac.utils import slice_channel, strip_invalid_lat, check_user_scanlines

LOG = logging.getLogger(__name__)

MISSING_DATA = -32001
MISSING_DATA_LATLON = -999999


def read_config():
    """Read output dir etc from config file."""
    conf = get_config()
    OUTDIR = conf.get('output', 'output_dir', raw=True)
    OUTPUT_FILE_PREFIX = conf.get('output', 'output_file_prefix', raw=True)

    SUNSATANGLES_DIR = os.environ.get('SM_SUNSATANGLES_DIR', OUTDIR)
    AVHRR_DIR = os.environ.get('SM_AVHRR_DIR', OUTDIR)
    QUAL_DIR = os.environ.get('SM_AVHRR_DIR', OUTDIR)

    return OUTPUT_FILE_PREFIX, SUNSATANGLES_DIR, AVHRR_DIR, QUAL_DIR


def save_gac(satellite_name,
             xutcs,
             lats, lons,
             ref1, ref2, ref3,
             bt3, bt4, bt5,
             sun_zen, sat_zen, sun_azi, sat_azi, rel_azi,
             qual_flags, start_line, end_line,
             gac_file, meta_data):
    corr = meta_data['sun_earth_distance_correction_factor']

    last_scan_line_number = qual_flags[-1, 0]

    # Strip invalid coordinates
    first_valid_lat, last_valid_lat = strip_invalid_lat(lats)
    if first_valid_lat > start_line:
        LOG.info('New start_line chosen (due to invalid lat/lon '
                 'info) = ' + str(first_valid_lat))
    if end_line > last_valid_lat:
        LOG.info('New end_line chosen (due to invalid lat/lon '
                 'info) = ' + str(last_valid_lat))

    # Check user-defined scanlines
    start_line, end_line = check_user_scanlines(
        start_line=start_line,
        end_line=end_line,
        first_valid_lat=first_valid_lat,
        last_valid_lat=last_valid_lat)

    # Slice data using new start/end lines
    ref1 = slice_channel(ref1,
                         start_line=start_line,
                         end_line=end_line,
                         first_valid_lat=first_valid_lat,
                         last_valid_lat=last_valid_lat)
    ref2 = slice_channel(ref2,
                         start_line=start_line,
                         end_line=end_line,
                         first_valid_lat=first_valid_lat,
                         last_valid_lat=last_valid_lat)
    ref3 = slice_channel(ref3,
                         start_line=start_line,
                         end_line=end_line,
                         first_valid_lat=first_valid_lat,
                         last_valid_lat=last_valid_lat)
    bt3 = slice_channel(bt3,
                        start_line=start_line,
                        end_line=end_line,
                        first_valid_lat=first_valid_lat,
                        last_valid_lat=last_valid_lat)
    bt4 = slice_channel(bt4,
                        start_line=start_line,
                        end_line=end_line,
                        first_valid_lat=first_valid_lat,
                        last_valid_lat=last_valid_lat)
    bt5 = slice_channel(bt5,
                        start_line=start_line,
                        end_line=end_line,
                        first_valid_lat=first_valid_lat,
                        last_valid_lat=last_valid_lat)
    sun_zen = slice_channel(sun_zen,
                            start_line=start_line,
                            end_line=end_line,
                            first_valid_lat=first_valid_lat,
                            last_valid_lat=last_valid_lat)
    sun_azi = slice_channel(sun_azi,
                            start_line=start_line,
                            end_line=end_line,
                            first_valid_lat=first_valid_lat,
                            last_valid_lat=last_valid_lat)
    sat_zen = slice_channel(sat_zen,
                            start_line=start_line,
                            end_line=end_line,
                            first_valid_lat=first_valid_lat,
                            last_valid_lat=last_valid_lat)
    sat_azi = slice_channel(sat_azi,
                            start_line=start_line,
                            end_line=end_line,
                            first_valid_lat=first_valid_lat,
                            last_valid_lat=last_valid_lat)
    rel_azi = slice_channel(rel_azi,
                            start_line=start_line,
                            end_line=end_line,
                            first_valid_lat=first_valid_lat,
                            last_valid_lat=last_valid_lat)
    lons = slice_channel(lons,
                         start_line=start_line,
                         end_line=end_line,
                         first_valid_lat=first_valid_lat,
                         last_valid_lat=last_valid_lat)
    lats = slice_channel(lats,
                         start_line=start_line,
                         end_line=end_line,
                         first_valid_lat=first_valid_lat,
                         last_valid_lat=last_valid_lat)
    qual_flags = slice_channel(qual_flags,
                               start_line=start_line,
                               end_line=end_line,
                               first_valid_lat=first_valid_lat,
                               last_valid_lat=last_valid_lat)
    xutcs = slice_channel(xutcs,
                          start_line=start_line,
                          end_line=end_line,
                          first_valid_lat=first_valid_lat,
                          last_valid_lat=last_valid_lat)

    total_number_of_scan_lines = lats.shape[0]

    # Reading time from the body of the gac file
    start = xutcs[0].astype(datetime.datetime)
    end = xutcs[-1].astype(datetime.datetime)
    startdate = start.strftime("%Y%m%d")
    starttime = start.strftime("%H%M%S%f")[:-5]
    enddate = end.strftime("%Y%m%d")
    endtime = end.strftime("%H%M%S%f")[:-5]

    # Apply scaling & offset
    bt3 -= 273.15
    bt4 -= 273.15
    bt5 -= 273.15
    for array in [bt3, bt4, bt5, ref1, ref2, ref3, sun_zen, sat_zen, sun_azi,
                  sat_azi, rel_azi]:
        array *= 100.0
    for array in [lats, lons]:
        array *= 1000.0

    # Replace NaN with fill values
    for array in [ref1, ref2, ref3, bt3, bt4, bt5, sun_zen, sat_zen, sun_azi,
                  sat_azi, rel_azi]:
        array[np.isnan(array)] = MISSING_DATA
    for array in [lats, lons]:
        array[np.isnan(array)] = MISSING_DATA_LATLON

    avhrrGAC_io(satellite_name, xutcs, startdate, enddate, starttime, endtime,
                lats, lons, ref1, ref2, ref3, bt3, bt4, bt5,
                sun_zen, sat_zen, sun_azi, sat_azi, rel_azi, qual_flags,
                start_line, end_line, total_number_of_scan_lines,
                last_scan_line_number, corr, gac_file)


def avhrrGAC_io(satellite_name, xutcs, startdate, enddate, starttime, endtime,
                arrLat_full, arrLon_full, ref1, ref2, ref3, bt3, bt4, bt5,
                arrSZA, arrSTZ, arrSAA, arrSTA, arrRAA, qual_flags,
                start_line, end_line, total_number_of_scan_lines,
                last_scan_line_number, corr, gac_file):
    import os

    # Read output dir etc from config file
    OUTPUT_FILE_PREFIX, SUNSATANGLES_DIR, AVHRR_DIR, QUAL_DIR = read_config()

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
    dset3 = g2.create_dataset("scanline_timestamps", dtype='int64',
                              data=xutcs.astype('int64'))
    dset3.attrs['units'] = 'Milliseconds since 1970-01-01 00:00:00 UTC'
    dset3.attrs['calendar'] = 'standard'

    fout.close()
