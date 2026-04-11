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

import csv
import unittest
from importlib.resources import files

import cftime
import numpy as np
import xarray as xr

from pygac.calibration.ir_uncertainty import (
    convBT,
    find_solar,
    get_gainval,
    get_uICT,
    get_uncert_parameter_thresholds,
    get_vars,
    ir_uncertainty,
    open_zenodo_uncert_file,
)
from pygac.calibration.noaa import Calibrator


#
# Code to read csv AVHRR data used for tests
#
def read_csv(filename):

    #
    # Read in data
    #
    d = np.loadtxt(filename,delimiter=",",skiprows=1)

    #
    # Mask
    #
    mask = np.zeros((d.shape[0]),dtype=np.int8)
    mask[:] = 0
    gd = (d[:,5] == 1)
    mask[gd] = 1

    bad_data = np.zeros((d.shape[0]),dtype=np.int8)
    bad_data[:] = 0
    gd = (d[:,6] == 1)
    bad_data[gd] = 1

    #
    # Get time units from csv header
    #
    with open(filename,"r") as fp:
        csv_reader = csv.reader(fp,delimiter=",")
        header = next(csv_reader)
    units = header[7]
    #
    # Get time to search stored values if needed
    # Convert to np.datetime64
    #
    times = cftime.num2date(d[:,7],units,
                            only_use_cftime_datetimes=False,
                            only_use_python_datetimes=True)
    outtime = []
    for i in range(len(times)):
        outtime.append(np.datetime64(times[i].strftime("%Y-%m-%dT%H:%M:%S.%f")))
    outtime = np.array(outtime)

    #
    # Create xarray data in form expected
    #
    da1 = xr.DataArray(name="scan_line_index",dims=["scan_line_index"],
                       data=d[:,0].astype(dtype=np.uint16))
    da2 = xr.DataArray(name="prt_counts",dims=["scan_line_index"],
                       data=d[:,1],attrs={"_FillValue":np.nan})
    data1 = np.zeros((d.shape[0],3))
    data1[:,0] = d[:,2]
    data1[:,1] = d[:,9]
    data1[:,2] = d[:,10]
    da3 = xr.DataArray(name="ict_counts",dims=["scan_line_index",
                                               "ir_channel_name"],
                       attrs={"_FillValue":np.nan},data=data1)
    data2 = np.zeros((d.shape[0],3))
    data2[:,0] = d[:,3]
    data2[:,1] = d[:,11]
    data2[:,2] = d[:,12]
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
#    da8 = xr.DataArray(name="times",dims=["scan_line_index"],
#                       attrs={"_FillValue":np.nan,
#                              "units":units},data=d[:,6])
    da8 = xr.DataArray(name="times",dims=["scan_line_index"],
                       data=outtime)
    da9 = xr.DataArray(name="uICT",dims=["scan_line_index"],
                       attrs={"_FillValue":np.nan},data=d[:,8])
    total_space = np.zeros((d.shape[0],10,6))
    total_space[:,:,0] = d[:,13:23]
    total_space[:,:,1] = d[:,23:33]
    total_space[:,:,2] = d[:,33:43]
    da10 = xr.DataArray(name="total_space_counts",
                        dims=["scan_line_index","ten_spots","channel_name"],
                        attrs={"_FillValue":np.nan},data=total_space)
    total_ict = np.zeros((d.shape[0],10,6))
    total_ict[:,:,0] = d[:,43:53]
    total_ict[:,:,1] = d[:,53:63]
    total_ict[:,:,2] = d[:,63:73]
    da11 = xr.DataArray(name="total_ict_counts",
                        dims=["scan_line_index","ten_spots","channel_name"],
                        attrs={"_FillValue":np.nan},data=total_ict)
    ds = xr.Dataset(data_vars=dict(channels=da6,counts=da7,
                                   prt_counts=da2,ict_counts=da3,
                                   space_counts=da4,sun_zen=da5,
                                   scan_line_index=da1,times=da8,
                                   uict=da9,total_space_counts=da10,
                                   total_ict_counts=da11),
                    attrs={"spacecraft_name":"noaa16"})
    
    return ds,mask,bad_data

#
# Input data from real GAC orbit with two solar peaks detected.
# Data taken from a section of real data and checked against code
#
#
# Test routine to get ICT uncertainty, check gain and solar detection as
# well
#
class TestGetUict(unittest.TestCase):

    def test_zenodo_uncert_file(self):

        platform = "metopa"
        with open_zenodo_uncert_file(platform, decode_times=False) as d:
            solar_start_time_2 = d["gain2_solar_start"].values[:]
            solar_stop_time_2 = d["gain2_solar_stop"].values[:]
        test_solar_start_time_2 = 1164115927.0
        test_solar_stop_time_2 = 1164116239.0
        np.testing.assert_allclose(solar_start_time_2[0],test_solar_start_time_2)
        np.testing.assert_allclose(solar_stop_time_2[0],test_solar_stop_time_2)

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
        ds,mask,bad_data = read_csv(csv_file)

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
        # Make solar flag plus bad data
        #
        solar_flag = bad_data[:]
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

        times = ds["times"].values[:]
        gd = np.isfinite(times)
        time = times[gd][0]

        #
        # Get gain at min PRT stddev point
        #
        gain_37,gain_time = get_gainval(time,times,avhrr_name,ict1,ict2,ict3,
                                        ict4,CS_1,CICT_1,CE_1,0.,mask,
                                        convT1,window,calculate=True)
        #
        # Check gain
        #
        check_gain = 0.002730006
        np.testing.assert_allclose(gain_37,check_gain,atol=0.00001)
        #
        # Get ICT uncertainty
        #
        uICT,nfigure = get_uICT(gain_37,CS_1,CICT_1,Tict,0.,convT1,bad_data,
                                solar_flag,window,plt=None,nfigure=1)
        #
        # Check uICT
        #
        gd = np.isfinite(uICT)&np.isfinite(ds["uict"].values)
        np.testing.assert_allclose(uICT[gd],ds["uict"].values[gd],atol=0.0001)

    def test_all_uncertainties(self):
        """Run complete IR uncertainty suite with a small amount of channel
        data to be checked against"""

        #
        # Read in orbits worth of calibration  data
        #
        csv_file = files("pygac") / "tests/test_ir_uncertainty_avhrr_data.csv"
        ds,mask,bad_data = read_csv(csv_file)

        #
        # Add in some channel counts
        #
        # Force all counts to a value of 800 and then put in real values
        # at a location of 13400:13410, 200:210
        starty = 13400
        endy = 13410
        startx = 200
        endx = 210
        #
        # 3.7 micron channel
        #
        ch37_counts = np.array([[708., 697., 694., 686., 699., 725., 746., 720., 746., 677.],
                                [672., 693., 729., 672., 672., 694., 703., 727., 772., 701.],
                                [662., 683., 718., 712., 709., 686., 696., 704., 770., 774.],
                                [753., 690., 709., 721., 728., 725., 685., 685., 704., 735.],
                                [771., 735., 714., 710., 715., 760., 702., 669., 660., 671.],
                                [756., 806., 707., 702., 715., 777., 702., 649., 649., 648.],
                                [837., 791., 735., 714., 723., 775., 691., 653., 658., 672.],
                                [823., 784., 759., 754., 765., 777., 689., 647., 674., 706.],
                                [806., 758., 699., 758., 772., 782., 726., 662., 670., 742.],
                                [791., 766., 725., 752., 784., 768., 741., 708., 650., 692.]])
        #
        # 11 micron channel
        #
        ch11_counts = np.array([[372., 373., 403., 364., 373., 382., 415., 454., 404., 361.],
                                [378., 370., 388., 363., 359., 370., 376., 420., 448., 397.],
                                [367., 360., 379., 391., 388., 367., 369., 374., 409., 428.],
                                [427., 375., 374., 379., 386., 379., 364., 359., 384., 388.],
                                [482., 423., 374., 371., 378., 399., 367., 356., 376., 365.],
                                [398., 495., 389., 373., 376., 417., 370., 352., 397., 354.],
                                [545., 504., 419., 379., 384., 422., 364., 356., 384., 361.],
                                [518., 482., 427., 396., 424., 471., 372., 350., 373., 374.],
                                [472., 466., 431., 425., 468., 477., 397., 351., 358., 392.],
                                [475., 474., 386., 400., 469., 460., 418., 389., 355., 366.]])
        #
        # 12 micron channel
        #
        ch12_counts = np.array([[382., 383., 408., 378., 384., 390., 418., 470., 412., 376.],
                                [391., 382., 394., 379., 374., 383., 389., 431., 450., 407.],
                                [381., 372., 387., 399., 396., 380., 382., 384., 408., 428.],
                                [428., 388., 383., 387., 393., 387., 377., 372., 398., 394.],
                                [470., 422., 383., 381., 386., 402., 378., 370., 399., 382.],
                                [401., 479., 398., 384., 386., 416., 381., 370., 425., 373.],
                                [532., 487., 420., 388., 391., 419., 374., 373., 408., 375.],
                                [504., 473., 425., 399., 424., 460., 381., 367., 385., 382.],
                                [460., 459., 429., 423., 458., 463., 399., 365., 371., 394.],
                                [466., 465., 393., 403., 458., 451., 416., 397., 370., 375.]])
        #
        # Add counts to xarray
        #
        counts = np.zeros((len(mask),409,6))
        counts[:,3:6] = 800.
        counts[starty:endy,startx:endx,3] = ch37_counts[:,:]
        counts[starty:endy,startx:endx,4] = ch11_counts[:,:]
        counts[starty:endy,startx:endx,5] = ch12_counts[:,:]
        da = xr.DataArray(name="counts",dims=["scan_line_index",
                                              "columns",
                                              "channel_name"],
                          data=counts)
        ds["counts"] = da
        #
        # Setup test values
        #
        # Random variables
        #
        ch3b_rand = [[0.06488484,0.06287615,0.06235265,0.06100407,0.06323078,0.06829536,0.0731193,0.06725055,
                      0.0731193,0.05956437],
                     [0.05880003,0.06218293,0.06916112,0.05880003,0.05880003,0.06235523,0.06395658,0.06872663,
                      0.08029056,0.06359265],
                     [0.057333,0.06051837,0.066846,0.0656554,0.06507797,0.0610071,0.06270365,0.06414084,
                      0.07968344,0.08090882],
                     [0.07490472,0.0616704,0.06507521,0.06745705,0.06894058,0.0682958,0.06084065,0.06084065,
                      0.06413812,0.0705002],
                     [0.07998042,0.07049857,0.06604254,0.06526487,0.06624026,0.07678931,0.06377035,0.0583475,
                      0.05704585,0.05864612],
                     [0.07569954,0.09248224,0.06469596,0.06377114,0.0662411,0.08185159,0.06377114,0.05554347,
                      0.05554347,0.05541134],
                     [0.10795469,0.08662684,0.07050246,0.06604606,0.06787558,0.08122112,0.06184123,0.05608146,
                      0.05676841,0.05879971],
                     [0.10030491,0.08416451,0.07651714,0.07516984,0.07820761,0.08185451,0.06150408,0.05528211,
                      0.05910371,0.06451102],
                     [0.09249021,0.07624613,0.06323454,0.07624613,0.08029267,0.08349333,0.06851309,0.05733332,
                      0.05850053,0.07214635],
                     [0.08662845,0.07849931,0.06829843,0.07464709,0.08416688,0.07908661,0.07190523,0.06488757,
                      0.05567839,0.0620116]]
        ch3b_rand = np.array(ch3b_rand)
        ch4_rand = [[0.08107775,0.08114852,0.08336898,0.08051869,0.08114852,0.08179459,0.08431342,0.08763486,
                     0.08344639,0.08031224],
                    [0.0815081,0.08093949,0.08223737,0.08045235,0.08017823,0.08093949,0.08136474,0.08471993,
                     0.08709979,0.08291199],
                    [0.08073146,0.08024838,0.08158201,0.0824622,0.08223929,0.08073146,0.08087122,0.0812241,
                     0.08384165,0.08538079],
                    [0.0852953,0.08129311,0.08122193,0.08157984,0.08208957,0.08157984,0.08052112,0.08017799,
                     0.08194288,0.0822371],
                    [0.09029119,0.08496328,0.08122037,0.08100805,0.08150629,0.0830625,0.08072774,0.07997282,
                     0.08136292,0.08058875],
                    [0.08298595,0.09161164,0.08230937,0.08114915,0.08136267,0.08447483,0.08093743,0.07970371,
                     0.08290989,0.07983777],
                    [0.09730858,0.09256386,0.0846385,0.0815799,0.08194295,0.08488296,0.08052117,0.07997443,
                     0.08194295,0.08031472],
                    [0.09410614,0.0902951,0.08529754,0.08283806,0.0850493,0.08922233,0.08108231,0.07957426,
                     0.08115309,0.08122406],
                    [0.0893152,0.08874462,0.08562922,0.08512892,0.08893355,0.0897995,0.0829111,0.07963816
                     ,0.08010933,0.08253413],
                    [0.08960327,0.08950641,0.08208743,0.08313842,0.08902695,0.08818372,0.08455525,0.08230905,
                     0.07990477,0.08065759]]
        ch4_rand = np.array(ch4_rand)
        ch5_rand = [[0.10022906,0.10031326,0.10249901,0.09989465,0.10039769,0.10090941,0.10341929,0.10869096,
                     0.10286383,0.09972887],
                    [0.10100141,0.10023485,0.10126136,0.09998367,0.09956975,0.10031906,0.10082935,0.10466409,
                     0.1065665,0.10241445],
                    [0.10015242,0.09940732,0.1006598,0.10170121,0.10143745,0.10006869,0.10023638,0.10040503,
                     0.10250652,0.10437539],
                    [0.10437179,0.10074176,0.10031714,0.10065635,0.10117252,0.10065635,0.09981551,0.09940392,
                     0.10160955,0.10125942],
                    [0.10869303,0.10379704,0.10031494,0.10014678,0.10056899,0.10196157,0.09989633,0.09923869,
                     0.10169551,0.10023074],
                    [0.10187376,0.10969971,0.10160847,0.1004005,0.1005701,0.10323598,0.1001479,0.09923979,
                     0.10408316,0.09948468],
                    [0.11630537,0.11062442,0.10361432,0.10074517,0.10100289,0.10352047,0.09957119,0.0994891,
                     0.10250652,0.09965351],
                    [0.11266802,0.10903441,0.10409106,0.10170433,0.10399578,0.10761974,0.10015545,0.09900441,
                     0.10049277,0.10023942],
                    [0.10761638,0.10751003,0.10447185,0.10389753,0.10740402,0.10793751,0.10170114,0.09884051,
                     0.09932556,0.10126281],
                    [0.10825939,0.10815094,0.1011737,0.10205421,0.10740164,0.10666908,0.10323831,0.10152283,
                     0.09924199,0.0996513]]
        ch5_rand = np.array(ch5_rand)
        #
        # Sytematic components
        #
        ch3b_sys = [[0.26468206,0.25793362,0.25618459,0.25169872,0.25912077,0.27626417,0.29287497,0.27270049,
                     0.29287497,0.24694288],
                    [0.24439431,0.25557727,0.27918007,0.24439431,0.24439431,0.2561519,0.26151367,0.27769248,
                     0.31791112,0.2602918],
                    [0.2395847,0.25002544,0.27125082,0.26721399,0.26526286,0.25164319,0.25728896,0.26210596,
                     0.31573641,0.32005674],
                    [0.29904745,0.25388407,0.26529536,0.27337336,0.27844018,0.27623483,0.25113004,0.25113004,
                     0.26213789,0.28379385],
                    [0.31687903,0.28382509,0.26860135,0.26596954,0.26927176,0.30566897,0.26093479,0.24295515,
                     0.23871648,0.24393213],
                    [0.30184835,0.36141419,0.26404043,0.26092872,0.26926543,0.32347498,0.26092872,0.2338594,
                     0.2338594,0.23343519],
                    [0.4174578,0.34039327,0.28379166,0.26856999,0.27479068,0.32120261,0.25444403,0.23556261,
                     0.23778155,0.2443978],
                    [0.3895697,0.33158893,0.30462871,0.29992502,0.31055157,0.3234023,0.25329496,0.23296904,
                     0.24536826,0.26335444],
                    [0.36138887,0.3037139,0.25909298,0.3037139,0.31792189,0.32924325,0.27696438,0.23961005,
                     0.24341496,0.28945877],
                    [0.34044321,0.31165264,0.27626736,0.29817648,0.33168069,0.31371755,0.28866419,0.26468418,
                     0.23429039,0.25504047]]
        ch3b_sys = np.array(ch3b_sys)
        ch4_sys = [[0.30530973,0.3052862,0.30458133,0.30549808,0.3052862,0.30507452,0.30429995,0.30338738,
                    0.30455787,0.30556875],
                   [0.30516814,0.30535637,0.30493308,0.3055212,0.30561545,0.30535637,0.30521519,0.30418237,
                    0.30352719,0.30472171],
                   [0.3054276,0.30559249,0.30514521,0.30486317,0.30493365,0.3054276,0.30538051,0.30526283,
                    0.30444073,0.30399556],
                   [0.30401814,0.30523839,0.30526191,0.30514431,0.30497975,0.30514431,0.30549732,0.30561511,
                    0.30502676,0.30493276],
                   [0.3027332,0.30411197,0.30526208,0.30533268,0.30516799,0.30467462,0.30542684,0.30568599,
                    0.30521503,0.30547393],
                   [0.30469826,0.3024302,0.3049096,0.30528578,0.3052152,0.30425268,0.30535638,0.30578048,
                    0.30472173,0.30573331],
                   [0.30126504,0.30222004,0.3042054,0.30514421,0.30502666,0.30413513,0.30549722,0.30568572,
                    0.30502666,0.30556789],
                   [0.30189407,0.30273319,0.30401833,0.30474508,0.30408857,0.30298986,0.30530919,0.30582754,
                    0.30528566,0.30526213],
                   [0.30296529,0.30310533,0.30392336,0.3040638,0.30305864,0.30284862,0.30472018,0.30580242,
                    0.30563739,0.30483757],
                   [0.30289547,0.3029188,0.30497871,0.30464998,0.30303549,0.30324561,0.30422794,0.30490823,
                    0.3057083,0.30544914]]
        ch4_sys = np.array(ch4_sys)
        ch5_sys = [[0.30491314,0.30488847,0.30427345,0.30501184,0.30486381,0.30471596,0.30402827,0.3027606,
                    0.30417532,0.30506122],
                   [0.30469088,0.30491268,0.30461703,0.3049867,0.30511016,0.30488801,0.30474013,0.30370974,
                    0.30324627,0.30429752],
                   [0.30493756,0.3051598,0.30478962,0.30449425,0.30456803,0.30496223,0.30491289,0.30486356,
                    0.30427318,0.30378326],
                   [0.30378261,0.30476426,0.3048875,0.3047889,0.30464114,0.3047889,0.30503555,0.30515906,
                    0.30451814,0.30461653],
                   [0.30275987,0.30392956,0.30488768,0.30493701,0.30481372,0.30442,0.30501104,0.30520868,
                    0.30449373,0.30491234],
                   [0.30444487,0.30254191,0.30451862,0.30486333,0.30481403,0.30407678,0.30493732,0.30520899,
                    0.30385643,0.30513484],
                   [0.30126293,0.30234783,0.30397843,0.30476438,0.30469049,0.30400292,0.30510977,0.30513448,
                    0.30427259,0.30508507],
                   [0.30193725,0.30268716,0.3038563,0.30449392,0.30388076,0.30300281,0.30493723,0.30528312,
                    0.30483858,0.30491256],
                   [0.30300155,0.30302585,0.30375713,0.30390387,0.30305017,0.30292865,0.3044925,0.3053311,
                    0.3051827,0.30461547],
                   [0.30285635,0.30288063,0.30464068,0.3043948,0.30305073,0.30322104,0.30407586,0.30454227,
                    0.30520803,0.30508447]]
        ch5_sys = np.array(ch5_sys)
        #
        # ch_to_ch covariance data
        #
        ch3b_ch2ch_covar = [[101,105,106,108,104,96,90,98,90,111],
                            [112,106,95,112,112,106,103,96,81,104],
                            [115,109,98,100,101,108,105,103,82,80],
                            [87,107,101,97,95,96,108,108,103,93],
                            [81,93,100,101,99,85,103,113,115,112],
                            [86,69,102,103,99,79,103,118,118,119],
                            [58,75,93,100,97,80,107,117,116,112],
                            [63,77,85,87,83,79,107,119,111,102],
                            [69,86,104,86,81,78,96,115,113,91],
                            [75,83,96,88,77,82,91,101,118,106]]
        ch4_ch2ch_covar = [[83,82,81,83,82,82,80,78,81,83],
                           [82,83,82,83,83,83,82,80,78,81],
                           [83,83,82,81,82,83,83,82,80,79],
                           [79,82,82,82,82,82,83,83,82,82],
                           [76,80,82,83,82,81,83,83,82,83],
                           [81,76,82,82,82,80,83,84,81,83],
                           [72,75,80,82,82,80,83,83,82,83],
                           [74,76,79,81,80,77,83,84,82,82],
                           [77,77,79,80,77,77,81,84,83,81],
                           [77,77,82,81,77,78,80,82,83,83]]
        ch5_ch2ch_covar = [[82,82,80,82,81,81,80,76,80,82],
                           [81,82,81,82,82,82,81,79,78,80],
                           [82,82,81,81,81,82,82,81,80,79],
                           [79,81,82,81,81,81,82,82,81,81],
                           [76,79,82,82,81,80,82,82,81,82],
                           [81,76,81,81,81,80,82,82,79,82],
                           [72,75,79,81,81,79,82,82,80,82],
                           [74,76,79,81,79,77,82,82,81,82],
                           [77,77,79,79,77,77,81,83,82,81],
                           [77,77,81,80,77,78,80,81,82,82]]
        ch3b_ch2ch_covar = np.array(ch3b_ch2ch_covar)
        ch4_ch2ch_covar = np.array(ch4_ch2ch_covar)
        ch5_ch2ch_covar = np.array(ch5_ch2ch_covar)
        #
        # Flags
        #
        #0-14 1
        #4088 1
        #6831 1
        #7096-7674 2
        #8776-9760 2
        #9928 1
        #13655-13668 1
        flags = np.zeros((13669),dtype=np.int32)
        flags[0:15] = 1
        flags[4088] = 1
        flags[6831] = 1
        flags[7096:7675] = 2
        flags[8776:9761] = 2
        flags[9928] = 1
        flags[13655:13669] = 1
        #
        # Run IR uncertainty codes
        #
        irdata = ir_uncertainty(ds,mask)

        #
        # Now do tests
        #
        # Random
        #
        np.testing.assert_allclose(irdata["random"].values[starty:endy,startx:endx,0],ch3b_rand[:,:],atol=0.0001)
        np.testing.assert_allclose(irdata["random"].values[starty:endy,startx:endx,1],ch4_rand[:,:],atol=0.0001)
        np.testing.assert_allclose(irdata["random"].values[starty:endy,startx:endx,2],ch5_rand[:,:],atol=0.0001)
        #
        # Systematic
        #
        np.testing.assert_allclose(irdata["systematic"].values[starty:endy,startx:endx,0],ch3b_sys[:,:],atol=0.0001)
        np.testing.assert_allclose(irdata["systematic"].values[starty:endy,startx:endx,1],ch4_sys[:,:],atol=0.0001)
        np.testing.assert_allclose(irdata["systematic"].values[starty:endy,startx:endx,2],ch5_sys[:,:],atol=0.0001)
        #
        # ch2ch covar
        #
        np.testing.assert_allclose(irdata["chan_covar_ratio"].values[starty:endy,startx:endx,0],ch3b_ch2ch_covar[:,:])
        np.testing.assert_allclose(irdata["chan_covar_ratio"].values[starty:endy,startx:endx,1],ch4_ch2ch_covar[:,:])
        np.testing.assert_allclose(irdata["chan_covar_ratio"].values[starty:endy,startx:endx,2],ch5_ch2ch_covar[:,:])
        #
        # Flags
        #
        np.testing.assert_allclose(irdata["uncert_flags"].values[:],flags[:])
