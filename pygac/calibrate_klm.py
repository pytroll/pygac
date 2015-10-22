#!/usr/bin/env python

# -*- coding: utf-8 -*-

# Copyright (c) 2014 Abhay Devasthale

# Author(s):

#   Abhay Devasthale <abhay.devasthale@smhi.se>

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


import numpy as np
import sys
MISSING_DATA = -32001


# calibrates solar channels of AVHRR


def calibrate_solar(counts, year, jday, spacecraft_id, channel3_switch, corr, number_of_data_records):

    # setting up calibration coefficients

    if spacecraft_id == 4:
                     # /*noaa-15*/
        Cs1 = 500.0
        Cs2 = 500.0
        Cs3 = 500.0
        Cdark1 = 39.0
        Cdark2 = 40.0
        Cdark3 = 39.0
        al1 = 0.0605
        bl1 = 0.447
        cl1 = -0.060
        ah1 = 0.1815
        bh1 = 0.447
        ch1 = -0.060
        al2 = 0.0675
        bl2 = 0.035
        cl2 = 0.007
        ah2 = 0.2025
        bh2 = 0.035
        ch2 = 0.007
        al3 = 0.0275
        bl3 = 0.00
        cl3 = 0.00
        ah3 = 0.1846
        bh3 = 0.00
        ch3 = 0.00
        Ldate = 1998.3641
    elif spacecraft_id == 2:
        # /*noaa-16*/
        Cs1 = 498.96
        Cs2 = 500.17
        Cs3 = 499.43
        Cdark1 = 39.3
        Cdark2 = 38.9
        Cdark3 = 38.4
        al1 = 0.0560
        bl1 = 0.306
        cl1 = 0.025
        ah1 = 0.1680
        bh1 = 0.306
        ch1 = 0.025
        al2 = 0.0580
        bl2 = 0.586
        cl2 = 0.036
        ah2 = 0.1740
        bh2 = 0.586
        ch2 = 0.036
        al3 = 0.0288
        bl3 = -0.81
        cl3 = 0.00
        ah3 = 0.2013
        bh3 = -0.81
        ch3 = 0.00
        Ldate = 2000.7228
    elif spacecraft_id == 6:
        # /*noaa-17*/
        Cs1 = 501.12
        Cs2 = 500.73
        Cs3 = 501.37
        Cdark1 = 39.99
        Cdark2 = 39.09
        Cdark3 = 42.09
        al1 = 0.0575
        bl1 = 1.707
        cl1 = -0.151
        ah1 = 0.1725
        bh1 = 1.707
        ch1 = -0.151
        al2 = 0.0650
        bl2 = 3.117
        cl2 = -0.265
        ah2 = 0.1950
        bh2 = 3.117
        ch2 = -0.265
        al3 = 0.0308
        bl3 = 4.06
        cl3 = -0.37
        ah3 = 0.2153
        bh3 = 4.06
        ch3 = -0.37
        Ldate = 2002.47912
    elif spacecraft_id == 7:
        # /*noaa-18*/
        Cs1 = 500.54
        Cs2 = 500.40
        Cs3 = 500.56
        Cdark1 = 39.44
        Cdark2 = 39.40
        Cdark3 = 37.51
        al1 = 0.0555
        bl1 = 3.068
        cl1 = -0.443
        ah1 = 0.1665
        bh1 = 3.068
        ch1 = -0.443
        al2 = 0.0595
        bl2 = 4.541
        cl2 = -0.611
        ah2 = 0.1785
        bh2 = 4.541
        ch2 = -0.611
        al3 = 0.0262
        bl3 = 0.00
        cl3 = 0.00
        ah3 = 0.1849
        bh3 = 0.00
        ch3 = 0.00
        Ldate = 2005.18891
    elif spacecraft_id == 8:
        # /*noaa-19*/
        Cs1 = 496.43
        Cs2 = 500.37
        Cs3 = 496.11
        Cdark1 = 38.8
        Cdark2 = 39.00
        Cdark3 = 39.4
        al1 = 0.0560
        bl1 = -5.985
        cl1 = -8.687
        ah1 = 0.1680
        bh1 = -5.985
        ch1 = -8.687
        al2 = 0.0585
        bl2 = 2.263
        cl2 = 0.748
        ah2 = 0.1755
        bh2 = 2.263
        ch2 = 0.748
        al3 = 0.0272
        bl3 = 0.00
        cl3 = 0.00
        ah3 = 0.1880
        bh3 = 0.00
        ch3 = 0.00
        Ldate = 2009.096
    elif spacecraft_id == 12:
        # /*metop-02*/
        Cs1 = 501.0
        Cs2 = 500.00
        Cs3 = 502.00
        Cdark1 = 40.43
        Cdark2 = 39.75
        Cdark3 = 41.80
        al1 = 0.0555
        bl1 = 1.797
        cl1 = -0.352
        ah1 = 0.1665
        bh1 = 1.797
        ch1 = -0.352
        al2 = 0.0635
        bl2 = 2.149
        cl2 = -0.225
        ah2 = 0.1905
        bh2 = 2.149
        ch2 = -0.225
        al3 = 0.0310
        bl3 = 4.11
        cl3 = 0.00
        ah3 = 0.2170
        bh3 = 4.11
        ch3 = 0.00
        Ldate = 2006.77
        # metop-01 is missing here. But as this code is never called, I do
        # not add it for now. /Sara Hornquist 2015-01-07
    else:
        print "wrong satellite id - exit"
        sys.exit(0)

    print 'year, jday, spacecraft-id, launch date - ', year, jday, spacecraft_id, Ldate

    t = (year + jday / 365.0) - Ldate

    # channel 1
    stl = 0.0
    sth = 0.0
    cindex = 0
    raw_counts = 0
    raw_counts = counts[:, 0::5]
    stl = (al1 * (100.0 + bl1 * t + cl1 * t * t)) / 100.0
    sth = (ah1 * (100.0 + bh1 * t + ch1 * t * t)) / 100.0
    r1 = (raw_counts - Cdark1) * stl
    cindex = np.where(raw_counts > Cs1)
    r1[cindex] = (Cs1 - Cdark1) * stl + (raw_counts[cindex] - Cs1) * sth
    r1 = r1 * corr

    # channel 2
    stl = 0.0
    sth = 0.0
    cindex = 0
    raw_counts = 0
    raw_counts = counts[:, 1::5]
    stl = (al2 * (100.0 + bl2 * t + cl2 * t * t)) / 100.0
    sth = (ah2 * (100.0 + bh2 * t + ch2 * t * t)) / 100.0
    r2 = (raw_counts - Cdark2) * stl
    cindex = np.where(raw_counts > Cs2)
    r2[cindex] = (Cs2 - Cdark2) * stl + (raw_counts[cindex] - Cs2) * sth
    r2 = r2 * corr

    # channel 3a (named as channel 6 in the output HDF5 file)
    iswitch = 0
    iswitch = np.where(channel3_switch == 1)
    r3 = np.ones((number_of_data_records, raw_counts.shape[1])) * MISSING_DATA

    if np.size(iswitch) > 0:
        stl = 0.0
        sth = 0.0
        cindex = 0
        raw_counts = 0
        raw_counts = counts[:, 2::5]
        stl = (al3 * (100.0 + bl3 * t + cl3 * t * t)) / 100.0
        sth = (ah3 * (100.0 + bh3 * t + ch3 * t * t)) / 100.0
        r3 = (raw_counts - Cdark3) * stl
        cindex = np.where(raw_counts > Cs3)
        r3[cindex] = (Cs3 - Cdark3) * stl + (raw_counts[cindex] - Cs3) * sth
        r3 = r3 * corr
        iswitch = 0
        iswitch = np.where(channel3_switch == 0)
        r3[iswitch, :] = MISSING_DATA
        iswitch = 0
        iswitch = np.where(channel3_switch == 2)
        r3[iswitch, :] = MISSING_DATA
    else:
        r3[:, :] = MISSING_DATA

    return r1, r2, r3


# calibrates thermal channels of AVHRR


def calibrate_thermal(raw_counts, prt, ict, space, number_of_data_records, spacecraft_id, channel, line_numbers):

    if spacecraft_id == 4:
                     # /*noaa-15*/
        d10 = 276.60157
        d11 = 0.051045
        d12 = 1.36328E-06
        d13 = 0.0
        d14 = 0.0
        d20 = 276.62531
        d21 = 0.050909
        d22 = 1.47266E-06
        d23 = 0.0
        d24 = 0.0
        d30 = 276.67413
        d31 = 0.050907
        d32 = 1.47656E-06
        d33 = 0.0
        d34 = 0.0
        d40 = 276.59258
        d41 = 0.050966
        d42 = 1.47656E-06
        d43 = 0.0
        d44 = 0.0
        if channel == 3:
            cWavenumber = 2695.9743
            A = 1.624481
            B = 1.0 / 1.001989
            Ns = 0.0
            b0 = 0.0
            b1 = 0.0
            b2 = 0.0
        elif channel == 4:
            cWavenumber = 925.4075
            A = 0.338243
            B = 1.0 / 1.001283
            Ns = -4.50
            b0 = 4.76
            b1 = -0.0932
            b2 = 0.0004524
        elif channel == 5:
            cWavenumber = 839.8979
            A = 0.304856
            B = 1.0 / 1.000977
            Ns = -3.61
            b0 = 3.83
            b1 = -0.0659
            b2 = 0.0002811

    elif spacecraft_id == 2:
                     # /*noaa-16*/
        d10 = 276.355
        d11 = 5.562E-02
        d12 = -1.590E-05
        d13 = 2.486E-08
        d14 = -1.199E-11
        d20 = 276.142
        d21 = 5.605E-02
        d22 = -1.707E-05
        d23 = 2.595E-08
        d24 = -1.224E-11
        d30 = 275.996
        d31 = 5.486E-02
        d32 = -1.223E-05
        d33 = 1.862E-08
        d34 = -0.853E-11
        d40 = 276.132
        d41 = 5.494E-02
        d42 = -1.344E-05
        d43 = 2.112E-08
        d44 = -1.001E-11
        if channel == 3:
            cWavenumber = 2681.2540
            A = 1.6774586
            B = 1.0 / 1.0017316
            Ns = 0.0
            b0 = 0.0
            b1 = 0.0
            b2 = 0.0
        elif channel == 4:
            cWavenumber = 922.34790
            A = 0.55636216
            B = 1.0 / 1.0014921
            Ns = -2.467
            b0 = 2.96
            b1 = -0.05411
            b2 = 0.00024532
        elif channel == 5:
            cWavenumber = 834.61814
            A = 0.41430789
            B = 1.0 / 1.0012166
            Ns = -2.009
            b0 = 2.25
            b1 = -0.03665
            b2 = 0.00014854

    elif spacecraft_id == 6:
                     # /*noaa-17*/
        d10 = 276.628
        d11 = 0.05098
        d12 = 1.371E-06
        d13 = 0.0
        d14 = 0.0
        d20 = 276.538
        d21 = 0.05098
        d22 = 1.371E-06
        d23 = 0.0
        d24 = 0.0
        d30 = 276.761
        d31 = 0.05097
        d32 = 1.369E-06
        d33 = 0.0
        d34 = 0.0
        d40 = 276.660
        d41 = 0.05100
        d42 = 1.348E-06
        d43 = 0.0
        d44 = 0.0
        if channel == 3:
            cWavenumber = 2669.1414
            A = 1.70002941
            B = 1.0 / 1.0026724
            Ns = 0.0
            b0 = 0.0
            b1 = 0.0
            b2 = 0.0
        elif channel == 4:
            cWavenumber = 928.29959
            A = 0.56634758
            B = 1.0 / 1.0015205
            Ns = -8.55
            b0 = 8.22
            b1 = -0.15795
            b2 = 0.00075579
        elif channel == 5:
            cWavenumber = 840.20289
            A = 0.37264803
            B = 1.0 / 1.0010841
            Ns = -3.97
            b0 = 4.31
            b1 = -0.07318
            b2 = 0.00030976

    elif spacecraft_id == 7:
                     # /*noaa-18*/
        d10 = 276.601
        d11 = 0.05090
        d12 = 1.657E-06
        d13 = 0.0
        d14 = 0.0
        d20 = 276.683
        d21 = 0.05101
        d22 = 1.482E-06
        d23 = 0.0
        d24 = 0.0
        d30 = 276.565
        d31 = 0.05117
        d32 = 1.313E-06
        d33 = 0.0
        d34 = 0.0
        d40 = 276.615
        d41 = 0.05103
        d42 = 1.484E-06
        d43 = 0.0
        d44 = 0.0
        if channel == 3:
            cWavenumber = 2660.6468
            A = 1.7222650
            B = 1.0 / 1.0028633
            Ns = 0.0
            b0 = 0.0
            b1 = 0.0
            b2 = 0.0
        elif channel == 4:
            cWavenumber = 928.73452
            A = 0.54696239
            B = 1.0 / 1.0014581
            Ns = -5.53
            b0 = 5.82
            b1 = -0.11069
            b2 = 0.00052337
        elif channel == 5:
            cWavenumber = 834.08306
            A = 0.39938376
            B = 1.0 / 1.0011724
            Ns = -2.22
            b0 = 2.67
            b1 = -0.04360
            b2 = 0.00017715

    elif spacecraft_id == 8:
                     # /*noaa-19*/
        d10 = 276.6067
        d11 = 0.051111
        d12 = 1.405783E-06
        d13 = 0.0
        d14 = 0.0
        d20 = 276.6119
        d21 = 0.051090
        d22 = 1.496037E-06
        d23 = 0.0
        d24 = 0.0
        d30 = 276.6311
        d31 = 0.051033
        d32 = 1.496990E-06
        d33 = 0.0
        d34 = 0.0
        d40 = 276.6268
        d41 = 0.051058
        d42 = 1.493110E-06
        d43 = 0.0
        d44 = 0.0
        if channel == 3:
            cWavenumber = 2670.2425
            A = 1.6863857
            B = 1.0 / 1.0025955
            Ns = 0.0
            b0 = 0.0
            b1 = 0.0
            b2 = 0.0
        elif channel == 4:
            cWavenumber = 927.92374
            A = 0.39419031
            B = 1.0 / 1.0013299
            Ns = -5.49
            b0 = 5.70
            b1 = -0.11187
            b2 = 0.00054668
        elif channel == 5:
            cWavenumber = 831.28619
            A = 0.26364620
            B = 1.0 / 1.0009546
            Ns = -3.39
            b0 = 3.58
            b1 = -0.05991
            b2 = 0.00024985

    elif spacecraft_id == 12:
                     # /*metop-02*/
        d10 = 276.6194
        d11 = 0.050919
        d12 = 1.471e-06
        d13 = 0.0
        d14 = 0.0
        d20 = 276.6511
        d21 = 0.050892
        d22 = 1.489e-06
        d23 = 0.0
        d24 = 0.0
        d30 = 276.6597
        d31 = 0.050845
        d32 = 1.521e-06
        d33 = 0.0
        d34 = 0.0
        d40 = 276.3685
        d41 = 0.050992
        d42 = 1.482e-06
        d43 = 0.0
        d44 = 0.0
        if channel == 3:
            cWavenumber = 2687.0392
            A = 2.0653147
            B = 1.0 / 1.0034418
            Ns = 0.0
            b0 = 0.0
            b1 = 0.0
            b2 = 0.0
        elif channel == 4:
            cWavenumber = 927.27630
            A = 0.56503332
            B = 1.0 / 1.0015090
            Ns = -4.98
            b0 = 5.44
            b1 = -0.10152
            b2 = 0.00046964
        elif channel == 5:
            cWavenumber = 837.80762
            A = 0.38472766
            B = 1.0 / 1.0011264
            Ns = -3.40
            b0 = 3.84
            b1 = -0.06249
            b2 = 0.00025239

    else:
        print "wrong satellite id - exit"
        sys.exit(0)

    # adjustment and preparation for calculating four PRT temperatures

    izeros = np.where(prt <= 50)
    izeros = izeros[0]

    inonzeros = np.where(prt > 50)
    inonzeros = inonzeros[0]

    offset = 0

    for i, prt_val in enumerate(prt):
        if prt_val < 0:
            offset = i
        break

    iprt = (line_numbers - 1 - offset) % 5

    tprt = np.zeros((float(number_of_data_records)))
    iones = np.where(iprt == 1)
    itwos = np.where(iprt == 2)
    ithrees = np.where(iprt == 3)
    ifours = np.where(iprt == 4)
    prt = prt.astype(float)

    tprt[iones] = d10 + d11 * prt[iones] + d12 * prt[iones] * prt[iones] + d13 * prt[iones] * \
        prt[iones] * prt[iones] + d14 * prt[iones] * \
        prt[iones] * prt[iones] * prt[iones]
    tprt[itwos] = d20 + d21 * prt[itwos] + d22 * prt[itwos] * prt[itwos] + d23 * prt[itwos] * \
        prt[itwos] * prt[itwos] + d24 * prt[itwos] * \
        prt[itwos] * prt[itwos] * prt[itwos]
    tprt[ithrees] = d30 + d31 * prt[ithrees] + d32 * prt[ithrees] * prt[ithrees] + d33 * prt[ithrees] * \
        prt[ithrees] * prt[ithrees] + d34 * prt[ithrees] * \
        prt[ithrees] * prt[ithrees] * prt[ithrees]
    tprt[ifours] = d40 + d41 * prt[ifours] + d42 * prt[ifours] * prt[ifours] + d43 * prt[ifours] * \
        prt[ifours] * prt[ifours] + d44 * prt[ifours] * \
        prt[ifours] * prt[ifours] * prt[ifours]

    utprt = tprt
    izeros = np.where(iprt == 0)
    izeros = izeros[0]

    inonzeros = np.where(iprt > 0)
    inonzeros = inonzeros[0]
    temp_prt = np.interp(izeros, inonzeros, tprt[inonzeros])
    tprt[izeros] = temp_prt

    # convolving and smoothing PRT, ICT and SPACE values

    if number_of_data_records > 51:
        window = 51
        # note that the window size has to be an odd number and greater than 2
    else:
        window = 3

    weighting_function = np.ones((window)) / window
    tprt_convolved = np.convolve(tprt, weighting_function, 'same')
    ict_convolved = np.convolve(ict, weighting_function, 'same')
    space_convolved = np.convolve(space, weighting_function, 'same')
    tprt_convolved[0:(window - 1) / 2] = tprt_convolved[(window - 1) / 2]
    ict_convolved[0:(window - 1) / 2] = ict_convolved[(window - 1) / 2]
    space_convolved[0:(window - 1) / 2] = space_convolved[(window - 1) / 2]
    tprt_convolved[-(window - 1) / 2:] = tprt_convolved[-((window + 1) / 2)]
    ict_convolved[-(window - 1) / 2:] = ict_convolved[-((window + 1) / 2)]
    space_convolved[-(window - 1) / 2:] = space_convolved[-((window + 1) / 2)]

    new_tprt = np.transpose(np.tile(tprt_convolved, (raw_counts.shape[1], 1)))
    new_ict = np.transpose(np.tile(ict_convolved, (raw_counts.shape[1], 1)))
    new_space = np.transpose(
        np.tile(space_convolved, (raw_counts.shape[1], 1)))

    # calibrating thermal channel

    tBB = new_tprt
    tsBB = A + B * tBB
    nBB = (1.1910427 * 0.000010) * cWavenumber * cWavenumber * cWavenumber
    nBB = nBB / (np.exp((1.4387752 * cWavenumber) / tsBB) - 1.0)

    Nlin = Ns + \
        (((nBB - Ns) * (new_space - raw_counts)) / (new_space - new_ict))
    Ncor = b0 + b1 * Nlin + b2 * Nlin * Nlin
    Ne = Nlin + Ncor
    tsE = 1.4387752 * cWavenumber
    tsE = tsE / \
        np.log(
            1.0 + ((1.1910427 * 0.000010 * cWavenumber * cWavenumber * cWavenumber) / Ne))
    bt = (tsE - A) / B

    return bt
