#!/usr/bin/env python

# Copyright (c) 2014-2025 Pytroll Developers

# Author(s):

#   Jonathan Mittaz <j.mittaz@reading.ac.uk>

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

"""Test function for the calibration PRT numbering for NOAA calibration
"""

import sys
import unittest

import numpy as np
import pandas as pd
import xarray as xr
import csv
import cftime
from datetime import datetime

from pygac.calibration.noaa import Calibrator
from pygac.calibration.ir_uncertainty import find_solar,\
    get_uncert_parameter_thresholds, get_vars, convBT,\
    get_gainval,get_uICT

#
# Code to read csv AVHRR data used for tests
#
def read_csv(filename):

    #
    # Read in data
    #
    d = np.loadtxt(filename,delimiter=',',skiprows=1)

    #
    # Mask
    #
    mask = np.zeros((d.shape[0]),dtype=bool)
    mask[:] = False
    gd = (d[:,5] == 1)
    mask[gd] = True

    #
    # Get time units from csv header
    #
    with open(filename,'r') as fp:
        csv_reader = csv.reader(fp,delimiter=',')
        header = next(csv_reader)
    units = header[6]
    #time = cftime.num2date(d[:,6],units,
    #                       only_use_cftime_datetimes=False,
    #                       only_use_python_datetimes=True)
    #dt64 = []
    #for i in range(len(time)):
    #    dt64.append(np.datetime64(time[i]))
    #dt64 = np.array(dt64)
    
    #
    # Create xarray data in form expected
    #
    da1 = xr.DataArray(name="scan_line_index",dims=["scan_line_index"],
                       data=d[:,0].astype(dtype=np.uint16))
    da2 = xr.DataArray(name="prt_counts",dims=["scan_line_index"],
                       data=d[:,1],attrs={"_FillValue":np.nan})
    data1 = np.zeros((d.shape[0],3))
    data1[:,0] = d[:,2]
    da3 = xr.DataArray(name="ict_counts",dims=["scan_line_index",
                                               "ir_channel_name"],
                       attrs={"_FillValue":np.nan},data=data1)
    data2 = np.zeros((d.shape[0],3))
    data2[:,0] = d[:,3]
    da4 = xr.DataArray(name="space_counts",dims=["scan_line_index",
                                                 "ir_channel_name"],
                       attrs={"_FillValue":np.nan},data=data2)
    data3 = np.zeros((d.shape[0],409))
    for i in range(d.shape[0]):
        data3[i,:] = d[i,4]
    da5 = xr.DataArray(name="sun_zen",dims=["scan_line_index",
                                            "columns"],
                       attrs={"_FillValue":np.nan},data=data3)
    data4 = np.zeros((d.shape[0],409,6))
    da6 = xr.DataArray(name="channels",dims=["scan_line_index",
                                             "columns",
                                             "channel_name"],
                       attrs={"_FillValue":np.nan},data=data4)
    data5 = np.zeros((d.shape[0],409,6))
    da7 = xr.DataArray(name="counts",dims=["scan_line_index",
                                           "columns",
                                           "channel_name"],
                       attrs={"_FillValue":np.nan},data=data5)
    da8 = xr.DataArray(name="times",dims=["scan_line_index"],
                       attrs={"_FillValue":np.nan,
                              "units":units},data=d[:,6])
    da9 = xr.DataArray(name="uICT",dims=["scan_line_index"],
                       attrs={"_FillValue":np.nan},data=d[:,7])
    ds = xr.Dataset(data_vars=dict(channels=da6,counts=da7,
                                   prt_counts=da2,ict_counts=da3,
                                   space_counts=da4,sun_zen=da5,
                                   scan_line_index=da1,times=da8,
                                   uict=da9),
                    attrs={"spacecraft_name":"noaa16"})
    
    return ds,mask

#
# Input data from real GAC orbit with two solar peaks detected.
# Data taken from a section of real data and checked against code
#
#
# Test routine to get ICT uncertainty, check gain and solar detection as
# well
#
class TestGetUict(unittest.TestCase):
    
    def test_get_Uict(self):

        #
        # Thresholds and kernels
        #
        window, prt_bias, prt_sys, prt_threshold, ict_threshold, \
            space_threshold = get_uncert_parameter_thresholds()
        #
        # Data taken from real case stored as CSV file
        #
        # Read in AVHRR data CSV file
        #
        csv_file = files("pygac") / "tests/test_ir_uncertainty_avhrr_data.csv"
        
        ds,mask = read_csv(csv_file)

        #
        # Setup Radiance to Temperature etc.
        #
        avhrr_name = ds.attrs["spacecraft_name"]
        cal = Calibrator(avhrr_name)
        convT1 = convBT(cal,0)

        #
        # Get variables
        #
        CS_1,CICT_1,CE_1,Tict,ict1,ict2,ict3,ict4 = get_vars(ds,0,convT1,
                                                             window,
                                                             prt_threshold,
                                                             ict_threshold,
                                                             space_threshold,
                                                             True,
                                                             cal,
                                                             mask,
                                                             out_prt=True)
        #
        # Make solar contamination flags
        #
        pos1_1,pos1_2,pos1,solza1,\
        pos2_1,pos2_2,pos2,solza2 = \
            find_solar(ds,mask,out_time=False,outgain=False)
        #
        # Make solar flag
        #
        solar_flag = np.zeros(CE_1.shape[0],dtype=np.uint8)
        if pos1_1 >= 0 and pos1_2 >= 0:
            solar_flag[pos1_1:pos1_2+1] = 1
        if pos2_1 >= 0 and pos2_2 >= 0:
            solar_flag[pos2_1:pos2_2+1] = 1
        
        #
        # Values to check against for solar detection
        #
        check_pos1_1 = 7096
        check_pos1_2 = 7674
        check_pos2_1 = 8776
        check_pos2_2 = 9760

        #
        # Do checks on solar detection
        #
        np.testing.assert_allclose(pos1_1,check_pos1_1,atol=2)
        np.testing.assert_allclose(pos1_2,check_pos1_2,atol=2)
        np.testing.assert_allclose(pos2_1,check_pos2_1,atol=2)
        np.testing.assert_allclose(pos2_2,check_pos2_2,atol=2)
        
        #
        # Get time to search stored values if needed
        # Convert to np.datetime64
        #
        times = cftime.num2date(ds['times'].values,ds['times'].attrs['units'],
                                only_use_cftime_datetimes=False,
                                only_use_python_datetimes=True)
                                
        gd = np.isfinite(ds['times'].values)
        time = times[gd][0]
        strtime = time.strftime('%Y-%m-%dT%H:%M:%S.%f'.
                                format(time.year,time.month,time.day,
                                       time.hour,time.minute,time.second,
                                       time.microsecond))
        time = np.datetime64(strtime)

        #
        # Get gain at min PRT stddev point
        #
        gain_37 = get_gainval(time,avhrr_name,ict1,ict2,ict3,ict4,CS_1,CICT_1,
                              CE_1,0.,mask,convT1,window,calculate=True)
        #
        # Check gain
        #
        check_gain = 0.0027243531
        np.testing.assert_allclose(gain_37,check_gain,atol=0.00001)
        #
        # Get ICT uncertainty
        #
        uICT,nfigure = get_uICT(gain_37,CS_1,CICT_1,Tict,0.,convT1,mask,
                                solar_flag,window,plt=None,nfigure=1)
        #
        # Check uICT
        #
        gd = np.isfinite(uICT)&np.isfinite(ds['uict'].values)
        np.testing.assert_allclose(uICT[gd],ds['uict'].values[gd],atol=0.0001)
        
