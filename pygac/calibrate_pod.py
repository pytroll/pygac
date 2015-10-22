#!/usr/bin/env python

# -*- coding: utf-8 -*-

# Copyright (c) 2014 Abhay Devasthale, Martin Raspaud

# Author(s):

#   Abhay Devasthale <abhay.devasthale@smhi.se>
#   Martin Raspaud   <martin.raspaud@smhi.se>

# This work was done in the framework of ESA-CCI-Clouds phase I


# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY
# without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import numpy as np

MISSING_DATA = -32001


# calibrates solar channels of AVHRR


def calibrate_solar_pod(counts, year, jday, spacecraft_id, channel3_switch, corr, number_of_data_records):

    # setting up calibration coefficients

    if spacecraft_id == 4:
        #/*noaa-7*/
        Cdark1 = 36.0
        Cdark2 = 37.0
        Cdark3 = 39.0

        al1 = 0.111206
        bl1 = 8.85810
        cl1 = -1.11938

        ah1 = 0.111206
        bh1 = 8.85810
        ch1 = -1.11938

        al2 = 0.107664
        bl2 = 20.7843
        cl2 = -4.32167

        ah2 = 0.107664
        bh2 = 20.7843
        ch2 = -4.32167

        al3 = 0.0000
        bl3 = 0.00
        cl3 = 0.00

        ah3 = 0.0000
        bh3 = 0.00
        ch3 = 0.00

        Ldate = 1981.4764

    elif spacecraft_id == 7:
        #/*noaa-9*/
        Cdark1 = 38.0
        Cdark2 = 40.0
        Cdark3 = 38.0

        al1 = 0.109312
        bl1 = 5.66526
        cl1 = 0.271520

        ah1 = 0.109312
        bh1 = 5.66526
        ch1 = 0.271520

        al2 = 0.116010
        bl2 = 3.93889
        cl2 = -0.258047

        ah2 = 0.116010
        bh2 = 3.93889
        ch2 = -0.258047

        al3 = 0.0000
        bl3 = 0.00
        cl3 = 0.00

        ah3 = 0.0000
        bh3 = 0.00
        ch3 = 0.00

        Ldate = 1984.9480

    elif spacecraft_id == 8:
        #/*noaa-10*/
        Cdark1 = 39.44
        Cdark2 = 39.40
        Cdark3 = 37.51

        al1 = 0.114803
        bl1 = 3.43456
        cl1 = -0.384635

        ah1 = 0.114803
        bh1 = 3.43456
        ch1 = -0.384635

        al2 = 0.132682
        bl2 = 0.00926103
        cl2 = 0.137333

        ah2 = 0.132682
        bh2 = 0.00926103
        ch2 = 0.137333

        al3 = 0.0
        bl3 = 0.00
        cl3 = 0.0

        ah3 = 0.0
        bh3 = 0.00
        ch3 = 0.0

        Ldate = 1986.712

    elif spacecraft_id == 1:
        #/*noaa-11*/
        Cdark1 = 40.0
        Cdark2 = 40.0
        Cdark3 = 40.0

        al1 = 0.114786
        bl1 = -1.82974
        cl1 = 0.412551

        ah1 = 0.114786
        bh1 = -1.82974
        ch1 = 0.412551

        al2 = 0.115219
        bl2 = -1.39495
        cl2 = 0.236978

        ah2 = 0.115219
        bh2 = -1.39495
        ch2 = 0.236978

        al3 = 0.0000
        bl3 = 0.00
        cl3 = 0.00

        ah3 = 0.0000
        bh3 = 0.00
        ch3 = 0.00

        Ldate = 1988.7064

    elif spacecraft_id == 5:
        # /*noaa-12*/
        Cdark1 = 41.0
        Cdark2 = 40.0
        Cdark3 = 40.0

        al1 = 0.118059
        bl1 = 4.88848
        cl1 = -0.347769

        ah1 = 0.118059
        bh1 = 4.88848
        ch1 = -0.347769

        al2 = 0.135684
        bl2 = 4.67098
        cl2 = -0.399447

        ah2 = 0.135684
        bh2 = 4.67098
        ch2 = -0.399447

        al3 = 0.000
        bl3 = 0.00
        cl3 = 0.00

        ah3 = 0.000
        bh3 = 0.00
        ch3 = 0.00

        Ldate = 1991.3669

    elif spacecraft_id == 3:
        # /*noaa-14*/
        Cdark1 = 41.0
        Cdark2 = 41.0
        Cdark3 = 39.0

        al1 = 0.121408
        bl1 = 4.07313
        cl1 = -0.375215

        ah1 = 0.121408
        bh1 = 4.07313
        ch1 = -0.375215

        al2 = 0.144515
        bl2 = 0.885331
        cl2 = 0.233727

        ah2 = 0.144515
        bh2 = 0.885331
        ch2 = 0.233727

        al3 = 0.0000
        bl3 = 0.00
        cl3 = 0.00

        ah3 = 0.0000
        bh3 = 0.00
        ch3 = 0.00

        Ldate = 1994.9699

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

    #r1=np.zeros((number_of_data_records, 409))

    raw_counts = counts[:, 0::5]

    stl = (al1 * (100.0 + bl1 * t + cl1 * t * t)) / 100.0
    # sth=(ah1*(100.0+bh1*t+ch1*t*t))/100.0;
    r1 = (raw_counts - Cdark1) * stl
    # cindex=np.where(raw_counts>Cs1);
    # r1[cindex]=(Cs1-Cdark1)*stl+(raw_counts[cindex]-Cs1)*sth;
    r1 = r1 * corr

    # channel 2
    stl = 0.0
    sth = 0.0
    cindex = 0
    raw_counts = 0

    #r2=np.zeros((number_of_data_records, 409));
    raw_counts = counts[:, 1::5]
    stl = (al2 * (100.0 + bl2 * t + cl2 * t * t)) / 100.0
    # sth=(ah2*(100.0+bh2*t+ch2*t*t))/100.0;
    r2 = (raw_counts - Cdark2) * stl
    # cindex=np.where(raw_counts>Cs2);
    # r2[cindex]=(Cs2-Cdark2)*stl+(raw_counts[cindex]-Cs2)*sth;
    r2 = r2 * corr

    # channel 3a (named as channel 6 in the output HDF5 file)
    r3 = np.zeros((number_of_data_records, raw_counts.shape[1]))
    r3[:, :] = MISSING_DATA

    return r1, r2, r3

 # calibrates thermal channels of AVHRR


def calibrate_thermal_pod(raw_counts, prt, ict, space, number_of_data_records, spacecraft_id, channel, line_numbers):

    if spacecraft_id == 4:
        # /*noaa-7*/
        d10 = 277.099
        d11 = 5.048E-2
        d12 = 2.823E-6
        d20 = 276.734
        d21 = 5.069E-2
        d22 = 2.493E-6
        d30 = 276.876
        d31 = 5.148E-2
        d32 = 1.040E-6
        d40 = 276.160
        d41 = 5.128E-2
        d42 = 1.414E-6
        if(channel == 3):
            Ns = 0.0
            cWavenumber = 2684.5233
            A = 1.94882690
            B = 1.0 / 1.0029260
            b1 = 1.0
            b2 = 0.0
            b0 = 0.0
        elif(channel == 4):
            Ns = -5.16
            cWavenumber = 928.23757
            A = 0.52807997
            B = 1.0 / 1.0014039
            b1 = 0.89783
            b2 = 0.0004819
            b0 = 5.25
        elif(channel == 5):
            Ns = -4.28
            cWavenumber = 841.52137
            A = 0.40557027
            B = 1.0 / 1.0011789
            b1 = 0.93683
            b2 = 0.0002425
            b0 = 3.93

    elif spacecraft_id == 7:
             # /*noaa-9*/
        d10 = 277.018000
        d11 = 0.051280
        d12 = 0.000000E-00
        d20 = 276.750000
        d21 = 0.051280
        d22 = 0.000000E-00
        d30 = 276.862000
        d31 = 0.051280
        d32 = 0.000000E-00
        d40 = 276.546000
        d41 = 0.051280
        d42 = 0.000000E-00
        if(channel == 3):
            Ns = 0.0
            cWavenumber = 2690.0451
            A = 1.8832662
            B = 1.0 / 1.0028978
            b1 = 1.0
            b2 = 0.0
            b0 = 0.0
        elif(channel == 4):
            Ns = -5.530
            cWavenumber = 930.50230
            A = 0.5115335
            B = 1.0 / 1.0013570
            b1 = 0.88643
            b2 = 0.0006033
            b0 = 5.24
        elif(channel == 5):
            Ns = -3.06
            cWavenumber = 845.75000
            A = 0.3882150
            B = 1.0 / 1.0011210
            b1 = 0.95311
            b2 = 0.0002198
            b0 = 2.42

    elif spacecraft_id == 8:
             # /*noaa-10*/
        d10 = 276.659
        d11 = 0.051275
        d12 = 1.363e-06
        d20 = 276.659
        d21 = 0.051275
        d22 = 1.363e-06
        d30 = 276.659
        d31 = 0.051275
        d32 = 1.363e-06
        d40 = 276.659
        d41 = 0.051275
        d42 = 1.363e-06
        if(channel == 3):
            Ns = 0.0
            cWavenumber = 2672.6392
            A = 1.80168640
            B = 1.0 / 1.0026383
            b1 = 1.0
            b2 = 0.0
            b0 = 0.0
        elif(channel == 4):
            Ns = 0.0
            cWavenumber = 910.51930
            A = 0.45707890
            B = 1.0 / 1.0012338
            b1 = 1.0
            b2 = 0.0
            b0 = 0.0
        elif(channel == 5):
            Ns = 0.0
            cWavenumber = 910.51930
            A = 0.45707890
            B = 1.0 / 1.0012338
            b1 = 1.0
            b2 = 0.0
            b0 = 0.0

    elif spacecraft_id == 1:
             # /*noaa-11*/
        d10 = 276.597
        d11 = 0.051275
        d12 = 1.363e-06
        d20 = 276.597
        d21 = 0.051275
        d22 = 1.363e-06
        d30 = 276.597
        d31 = 0.051275
        d32 = 1.363e-06
        d40 = 276.597
        d41 = 0.051275
        d42 = 1.363e-06
        if(channel == 3):
            Ns = 0.0
            cWavenumber = 2680.05
            A = 1.738973
            B = 1.0 / 1.003354
            b1 = 1.0
            b2 = 0.0
            b0 = 0.0
        elif(channel == 4):
            Ns = -8.055
            cWavenumber = 927.462
            A = 0.321199
            B = 1.0 / 1.001213
            b1 = 0.84120
            b2 = 0.0008739
            b0 = 7.21
        elif(channel == 5):
            Ns = -3.51
            cWavenumber = 840.746
            A = 0.048652
            B = 1.0 / 1.000664
            b1 = 0.94598
            b2 = 0.0002504
            b0 = 2.92

    elif spacecraft_id == 5:
             # /*noaa-12*/
        d10 = 276.597
        d11 = 0.051275
        d12 = 1.363e-06
        d20 = 276.597
        d21 = 0.051275
        d22 = 1.363e-06
        d30 = 276.597
        d31 = 0.051275
        d32 = 1.363e-06
        d40 = 276.597
        d41 = 0.051275
        d42 = 1.363e-06
        if(channel == 3):
            Ns = 0.0
            cWavenumber = 2651.7708
            A = 1.90527390
            B = 1.0 / 1.0030100
            b1 = 1.0
            b2 = 0.0
            b0 = 0.0
        elif(channel == 4):
            Ns = -5.510
            cWavenumber = 922.36261
            A = 0.63404209
            B = 1.0 / 1.0017076
            b1 = 0.88929
            b2 = 0.0005968
            b0 = 5.11
        elif(channel == 5):
            Ns = -2.51
            cWavenumber = 838.02678
            A = 0.41086587
            B = 1.0 / 1.0012010
            b1 = 0.96299
            b2 = 0.0001775
            b0 = 1.91

    elif spacecraft_id == 3:
             # /*noaa-14*/
        d10 = 276.597
        d11 = 0.051275
        d12 = 1.363e-06
        d20 = 276.597
        d21 = 0.051275
        d22 = 1.363e-06
        d30 = 276.597
        d31 = 0.051275
        d32 = 1.363e-06
        d40 = 276.597
        d41 = 0.051275
        d42 = 1.363e-06
        if(channel == 3):
            Ns = 0.0069
            cWavenumber = 2654.25
            A = 1.885330
            B = 1.0 / 1.003839
            b1 = 1.00359
            b2 = 0.0
            b0 = -0.0031
        elif(channel == 4):
            Ns = -4.05
            cWavenumber = 928.349
            A = 0.308384
            B = 1.0 / 1.001443
            b1 = 0.92378
            b2 = 0.000382
            b0 = 3.72
        elif(channel == 5):
            Ns = -2.29
            cWavenumber = 833.040
            A = 0.022171
            B = 1.0 / 1.000538
            b1 = 0.96194
            b2 = 0.0001742
            b0 = 2.00
    else:
        print "wrong satellite id - exit"
        sys.exit(0)

    columns = raw_counts.shape[1]

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

    tprt = np.zeros((int(number_of_data_records)))
    iones = np.where(iprt == 1)
    itwos = np.where(iprt == 2)
    ithrees = np.where(iprt == 3)
    ifours = np.where(iprt == 4)
    prt = prt.astype(float)
    tprt[iones] = d10 + d11 * prt[iones] + d12 * (prt[iones] * prt[iones])
    tprt[itwos] = d20 + d21 * prt[itwos] + d22 * (prt[itwos] * prt[itwos])
    tprt[ithrees] = d30 + d31 * prt[ithrees] + \
        d32 * (prt[ithrees] * prt[ithrees])
    tprt[ifours] = d40 + d41 * prt[ifours] + d42 * (prt[ifours] * prt[ifours])

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

    weighting_function = np.ones(window, dtype=float) / window
    tprt_convolved = np.convolve(tprt, weighting_function, 'same')
    ict_convolved = np.convolve(ict, weighting_function, 'same')
    space_convolved = np.convolve(space, weighting_function, 'same')
    print tprt
    print tprt_convolved
    tprt_convolved[0:(window - 1) / 2] = tprt_convolved[(window - 1) / 2]
    ict_convolved[0:(window - 1) / 2] = ict_convolved[(window - 1) / 2]
    space_convolved[0:(window - 1) / 2] = space_convolved[(window - 1) / 2]
    tprt_convolved[-(window - 1) / 2:] = tprt_convolved[-((window + 1) / 2)]
    ict_convolved[-(window - 1) / 2:] = ict_convolved[-((window + 1) / 2)]
    space_convolved[-(window - 1) / 2:] = space_convolved[-((window + 1) / 2)]
    print tprt_convolved
    new_tprt = np.transpose(np.tile(tprt_convolved, (columns, 1)))
    new_ict = np.transpose(np.tile(ict_convolved, (columns, 1)))
    new_space = np.transpose(np.tile(space_convolved, (columns, 1)))

    # calibrating thermal channel

    tBB = new_tprt
    tsBB = A + B * tBB
    nBB = (1.1910427 * 0.000010) * cWavenumber * cWavenumber * cWavenumber
    nBB = nBB / (np.exp((1.4387752 * cWavenumber) / tsBB) - 1.0)

    Nlin = Ns + \
        (((nBB - Ns) * (new_space - raw_counts)) / (new_space - new_ict))
    Ncor = b0 + b1 * Nlin + b2 * Nlin * Nlin
    Ne = Ncor
    tsE = 1.4387752 * cWavenumber
    tsE = tsE / \
        np.log(
            1.0 + ((1.1910427 * 0.000010 * cWavenumber * cWavenumber * cWavenumber) / Ne))
    bt = (tsE - A) / B

    return bt
