#!/usr/bin/env python

# -*- coding: utf-8 -*-

# Copyright (c) 2014 Abhay Devasthale

# Author(s):

#   Abhay Devasthale <abhay.devasthale@smhi.se>
#   Adam Dybbroe <adam.dybbroe@smhi.se>

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

import sys
import numpy as np
import geotiepoints as gtp
import datetime
from pyorbital import astronomy
import pygac.calibrate_pod as cal_pod
from pygac import gac_io
import pprint

from pyorbital.orbital import Orbital

import logging
LOG = logging.getLogger(__name__)

import ConfigParser
import os

try:
    CONFIG_FILE = os.environ['PYGAC_CONFIG_FILE']
except KeyError:
    LOG.exception('Environment variable PYGAC_CONFIG_FILE not set!')
    raise

if not os.path.exists(CONFIG_FILE) or not os.path.isfile(CONFIG_FILE):
    raise IOError(str(CONFIG_FILE) + " pointed to by the environment " +
                  "variable PYGAC_CONFIG_FILE is not a file or does not exist!")


# setting up constants used in the processing

AVHRR_SWATH_WIDTH_DEG=55.385
sat_altitude=850000.0
M_PI=3.14159265359
AHA_DEG_TO_RAD=M_PI/180.
AHA_RAD_TO_DEG=180./M_PI
earth_radius=6370997.0
wholeday_ms=24.0*3600.0*1000.0
GAC_NEAR_NADIR_POSITION=207
MISSING_DATA = -32001


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
                   ("reserved_for_mounting_and_fixed_attitude_correction_indicator", ">i1"),
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
                     ("time_code", ">u1", (6, )),
                     ("quality_indicators", ">u4"),
                     ("calibration_coefficients", ">i4", (10, )),
                     ("number_of_meaningful_zenith_angles_and_earth_location_appended", ">u1"),
                     ("solar_zenith_angles", "i1", (51, )),
                     ("earth_location", ">i2", (102, )),
                     ("telemetry", ">u4", (35, )),
                     ("sensor_data", ">u4", (682, )),
                     ("add_on_zenith", ">u2", (10, )),
                     ("clock_drift_delta", ">u2"),
                     ("spare3", "u2", (11, ))])


def main(filename):

    with open(filename) as fd_:
        head = np.fromfile(fd_, dtype=header, count=1)
        scans = np.fromfile(fd_, dtype=scanline, count=head["number_of_scans"])

    # cleaning up the data
    if scans["scan_line_number"][0] != 1:
        scans = scans[1:]
    scans = scans[scans["scan_line_number"] != 0]

    number_of_scans = len(scans)

    spacecraft_id=int(head["noaa_spacecraft_identification_code"])
    spacecrafts = {4: 'noaa07',
                   7: 'noaa09',
		   8: 'noaa10',
		   1: 'noaa11',
		   5: 'noaa12',
		   3: 'noaa14',
                   }

    try:
        satellite_name = spacecrafts[spacecraft_id]
    except KeyError:
        raise KeyError("Wrong satellite id: " + str(spacecraft_id))

    conf = ConfigParser.ConfigParser()
    try:
        conf.read(CONFIG_FILE)
    except ConfigParser.NoSectionError:
        LOG.exception('Failed reading configuration file: ' + str(CONFIG_FILE))
        raise

    values = {"satname": satellite_name, }

    options = {}
    for option, value in conf.items('tle', raw = True):
        options[option] = value

    tle_filename = os.path.join(options['tledir'],
                                options["tlename"] % values)
    LOG.info('TLE filename = ' + str(tle_filename))


    # unpacking raw data to get counts

    packed_data = scans["sensor_data"]
    LOG.debug("Header: %s", pprint.pformat(dict(zip(head[0].dtype.names,
                                                    head[0]))))
    gac_counts = np.zeros((int(number_of_scans), 409 * 5))
    gac_counts[:, 0::3] = (packed_data & (1023 << 20)) >> 20
    gac_counts[:, 1::3] = (packed_data & (1023 << 10)) >> 10
    gac_counts[:, 2::3] = (packed_data & 1023)[:, :-1]



    # interpolating lat-on points using PYTROLL geotiepoints

    arrLat = np.zeros((int(number_of_scans), 51))
    arrLon = np.zeros((int(number_of_scans), 51))
    arrLat = scans["earth_location"][:, 0::2] / 128.0
    arrLon = scans["earth_location"][:, 1::2] / 128.0


    arrLon_full, arrLat_full = gtp.Gac_Lat_Lon_Interpolator(arrLon,arrLat)


    # getting time information and calculating solar zenith angle for the entire orbit

    arrYear = ((np.uint16(scans["time_code"][:,0]) << 8) | (np.uint16(scans["time_code"][:,1]) & 65535)) >> 9
    arrYear = np.where(arrYear>75, arrYear+1900, arrYear+2000)
    arrJDay = ((np.uint16(scans["time_code"][:,0]) << 8) | (np.uint16(scans["time_code"][:,1]) & 65535)) & 0x01ff

    arrUTC = ((np.uint32(scans["time_code"][:,2] & 7) << 24) | (np.uint32(scans["time_code"][:,3]) << 16) | (np.uint32(scans["time_code"][:,4]) << 8) | (np.uint32(scans["time_code"][:,5])));
    #arrUTC = ((np.uint32(scans["time_code"][:,2]) << 24) | (np.uint32(scans["time_code"][:,3]) << 16) | (np.uint32(scans["time_code"][:,4]) << 8) | (np.uint32(scans["time_code"][:,5]) & 4294967295)) & 0x07ffffff;
    #print arrYear, arrJDay, arrUTC

    # Calculating solar zenith angle

    arrSZA = np.zeros((int(number_of_scans), 409))
    for i in range(number_of_scans):
            temp_utc = datetime.datetime(int(arrYear[i]), 1, 1) + datetime.timedelta(int(arrJDay[i]) - 1) + datetime.timedelta(milliseconds=int(arrUTC[i]))
            arrSZA[i,:]= astronomy.sun_zenith_angle(temp_utc, arrLon_full[i,:], arrLat_full[i,:])



    # calculating satellite zenith angle

    pixel_pos = np.zeros((int(number_of_scans), 409))
    for i in range(number_of_scans):
            pixel_pos[i,:] = np.arange(0,409,1)

    scan_angle = AVHRR_SWATH_WIDTH_DEG*np.divide(np.absolute(pixel_pos-205),205.0)
    arrSTZ = np.zeros((int(number_of_scans), 409))
    arrSTZ = np.arcsin((1.0+sat_altitude/earth_radius)*np.sin(scan_angle*AHA_DEG_TO_RAD))*AHA_RAD_TO_DEG;



    # calculating solar azimuth angle

    arrSAA = np.zeros((int(number_of_scans), 409))
    for i in range(number_of_scans):
            temp_utc = datetime.datetime(int(arrYear[i]), 1, 1) + datetime.timedelta(int(arrJDay[i]) - 1) + datetime.timedelta(milliseconds=int(arrUTC[i]));
            temp_alt, temp_suna = astronomy.get_alt_az(temp_utc, arrLon_full[i,:], arrLat_full[i,:]);
            temp_suna=temp_suna*180.0/M_PI;
            suna=temp_suna*0.0;
            ii=np.where(temp_suna<0.0);
            suna[ii]=temp_suna[ii]+180.0;
            jj=np.where(temp_suna>=0.0);
            suna[jj]=temp_suna[jj]-180.0;
            arrSAA[i,:]=suna;


    # calculating satellite azimuth angle

    tle_name = tle_filename # './noaa_tle/TLE_%s.txt' % satellite_name
    if arrYear[1]<=1999:
	yearJday = '%02d%03d.' % (arrYear[1]%1900,arrJDay[1])
    else:
	yearJday = '%02d%03d.' % (arrYear[1]%2000,arrJDay[1])


    tlep=open(tle_name,'rt')
    tle_data=tlep.readlines()
    tlep.close()

    def find_tle_index(tle_data, yearJday):
    	for iindex, stringYearJday in enumerate(tle_data):
        	if yearJday in stringYearJday:
              		return iindex
    	return -1

    try:
    	iindex=find_tle_index(tle_data,yearJday)
    	tle1=tle_data[iindex]
    	tle2=tle_data[iindex+1]
    except KeyError:
    	print "Could't find TLE - exit";
        sys.exit(1);


    spacecraft_id=int(head["noaa_spacecraft_identification_code"])
    spacecrafts_orbital = {4: 'noaa 7',
                   7: 'noaa 9',
		   8: 'noaa 10',
		   1: 'noaa 11',
		   5: 'noaa 12',
		   3: 'noaa 14',
                   }

    try:
        satellite_name_orbital = spacecrafts_orbital[spacecraft_id]
    except KeyError:
        print "wrong satellite id - exit";
        sys.exit(1);


    orb = Orbital(satellite_name_orbital, line1=tle1, line2=tle2)

    arrSTA = np.zeros((int(number_of_scans), 409), dtype='float')
    for i in range(number_of_scans):
            temp_utc = datetime.datetime(int(arrYear[i]), 1, 1) + datetime.timedelta(int(arrJDay[i]) - 1) + datetime.timedelta(milliseconds=int(arrUTC[i]));
    	    arr_psta, arr_ele = orb.get_observer_look(temp_utc, arrLon_full[i,:], arrLat_full[i,:], 0)
    	    arrSTA[i,:]=arr_psta;

    arrSTA[:,205]=arrSTA[:,204]



    # calculating relative azimuth angle

    arrRAA = np.zeros((int(number_of_scans), 409));
    arrRAA = np.absolute(arrSTA-arrSAA);
    arrRAA = np.where(arrRAA>180.0, 360.0-arrRAA, arrRAA);


    # scaling angles and lat-lon values


    arrSTA=arrSTA-180.0;
    arrRAA=np.where(arrRAA<0.0,-1.0*arrRAA,arrRAA);
    arrRAA=180.0-arrRAA;

    arrSZA=arrSZA*100.0;
    arrSTZ=arrSTZ*100.0;
    arrSAA=arrSAA*100.0;
    arrSTA=arrSTA*100.0;
    arrRAA=arrRAA*100.0;

    arrLat_full=arrLat_full*100.0;
    arrLon_full=arrLon_full*100.0;


    # Earth-Sun distance correction factor

    corr = 1.0 - 0.0334*np.cos(2.0*M_PI*(arrJDay[0]-2)/365.25)



    # Calibrating solar channels

    channel3_switch=np.zeros(int(number_of_scans));
    channel3_switch=0;
    ref1,ref2,ref3=cal_pod.calibrate_solar_pod(gac_counts, int(arrYear[0]), int(arrJDay[0]), int(head["noaa_spacecraft_identification_code"]), channel3_switch, corr, int(number_of_scans));


    # Calibrating thermal channels

    decode_tele=np.zeros((int(number_of_scans),140));
    j=0;
    for i in range (35):
            decode_tele[:,j] = (scans["telemetry"][:,i] << 2) >> 22;
            j=j+1;
            decode_tele[:,j] = (scans["telemetry"][:,i] << 12) >> 22;
            j=j+1;
            decode_tele[:,j] = (scans["telemetry"][:,i] << 22) >> 22;
            j=j+1;


    prt_counts=np.mean(decode_tele[:,17:20], axis=1);


    # getting ICT counts

    ict_counts=np.zeros((int(number_of_scans),3));
    ict_counts[:,0]=np.mean(decode_tele[:,22:50:3],axis=1);
    ict_counts[:,1]=np.mean(decode_tele[:,23:51:3],axis=1);
    ict_counts[:,2]=np.mean(decode_tele[:,24:52:3],axis=1);


    # getting space counts

    space_counts=np.zeros((int(number_of_scans),3));
    space_counts[:,0]=np.mean(decode_tele[:,54:100:5],axis=1);
    space_counts[:,1]=np.mean(decode_tele[:,55:101:5],axis=1);
    space_counts[:,2]=np.mean(decode_tele[:,56:102:5],axis=1);



    # calibrating channels 3b, 4 and 5

    bt3=cal_pod.calibrate_thermal_pod(gac_counts[:,2::5],
                                      prt_counts,
                                      ict_counts[:,0],
                                      space_counts[:,0],
                                      int(number_of_scans),
                                      int(head["noaa_spacecraft_identification_code"]),
                                      channel=3,
                                      line_numbers=scans["scan_line_number"])
    bt4=cal_pod.calibrate_thermal_pod(gac_counts[:,3::5],
                                      prt_counts,
                                      ict_counts[:,1],
                                      space_counts[:,1],
                                      int(number_of_scans),
                                      int(head["noaa_spacecraft_identification_code"]),
                                      channel=4,
                                      line_numbers=scans["scan_line_number"])
    bt5=cal_pod.calibrate_thermal_pod(gac_counts[:,4::5],
                                      prt_counts,
                                      ict_counts[:,2],
                                      space_counts[:,2],
                                      int(number_of_scans),
                                      int(head["noaa_spacecraft_identification_code"]),
                                      channel=5,
                                      line_numbers=scans["scan_line_number"])


    bt3=(bt3-273.15)*100.0;
    bt4=(bt4-273.15)*100.0;
    bt5=(bt5-273.15)*100.0;


    # masking out corrupt scanlines

    scanline_quality=np.zeros((int(number_of_scans),3));
    scanline_quality[:,0]=scans["quality_indicators"]>>31;
    scanline_quality[:,1]=(scans["quality_indicators"]<<4)>>31;
    scanline_quality[:,2]=(scans["quality_indicators"]<<5)>>31;

    ii=np.where((scanline_quality[:,0]==1) | (scanline_quality[:,1]==1) | (scanline_quality[:,2]==1));
    arrLat_full[ii]=MISSING_DATA;
    arrLon_full[ii]=MISSING_DATA;
    ref1[ii]=MISSING_DATA;
    ref2[ii]=MISSING_DATA;
    ref3[ii]=MISSING_DATA;
    bt3[ii]=MISSING_DATA;
    bt4[ii]=MISSING_DATA;
    bt5[ii]=MISSING_DATA;
    arrSZA[ii]=MISSING_DATA;
    arrSTZ[ii]=MISSING_DATA;
    arrSAA[ii]=MISSING_DATA;
    arrSTA[ii]=MISSING_DATA;
    arrRAA[ii]=MISSING_DATA;


    # writing out calibrated AVHRR channel data and various sun-sat angles


    t = datetime.datetime(int(arrYear[1]), 1, 1) + datetime.timedelta(int(arrJDay[1]) - 1) + datetime.timedelta(milliseconds=int(arrUTC[1]))
    tenth_s = int(t.microsecond/100000)
    startdate = '%d%02d%02d' % (t.year,t.month,t.day)
    starttime = '%02d%02d%02d%01d' % (t.hour,t.minute,t.second,tenth_s)

    t = datetime.datetime(int(arrYear[-1]), 1, 1) + datetime.timedelta(int(arrJDay[-1]) - 1) + datetime.timedelta(milliseconds=int(arrUTC[-1]))
    tenth_s = int(t.microsecond/100000)
    enddate = '%d%02d%02d' % (t.year,t.month,t.day)
    endtime = '%02d%02d%02d%01d' % (t.hour,t.minute,t.second,tenth_s)

    gac_io.avhrrGAC_io(satellite_name, startdate, enddate, starttime, endtime, arrLat_full, arrLon_full, ref1, ref2, ref3, bt3, bt4, bt5, arrSZA, arrSTZ, arrSAA, arrSTA, arrRAA)



if __name__ == "__main__":
    main(filename)


