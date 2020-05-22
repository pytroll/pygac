#!/usr/bin/env python

# Copyright (c) 2014-2015, 2019 Pytroll Developers

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
#   Abhay Devasthale <abhay.devasthale@smhi.se>
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

"""Calibration coefficients and generic calibration functions
"""
from __future__ import division
import logging
import numpy as np
import json
from collections import namedtuple
from pkg_resources import resource_filename

import pygac.configuration

LOG = logging.getLogger(__name__)

class Calibrator(object):
    """Factory class to create namedtuples holding the calibration coefficients.

    Attributes:
        fields: coefficient names
        Calibrator: namedtuple constructor
        default_coeffs: dictonary containing default values for all spacecrafts
    """
    fields = 'ah al bh bl ch cl c_s c_dark l_date d n_s c_wn a b b0 b1 b2'.split()
    Calibrator = namedtuple('Calibrator', fields)
    default_coeffs = None

    def __new__(cls, spacecraft, custom_coeffs=None):
        """Creates a namedtuple for calibration coefficients of a given spacecraft

        Args:
            spacecraft (str): spacecraft name in pygac convention
            custom_coeffs (dict): custom coefficients (optional)

        Returns:
            calibrator (namedtuple): calibration coefficients
        """
        if cls.default_coeffs is None:
            cls.load_defaults()
        custom_coeffs = custom_coeffs or {}
        customs = {key: cls.parse(value) for key, value in custom_coeffs.items()}
        defaults = cls.default_coeffs[spacecraft]
        spacecraft_coeffs = dict.fromkeys(cls.fields)
        spacecraft_coeffs.update(defaults)
        spacecraft_coeffs.update(customs)
        if custom_coeffs:
            LOG.info('Using following custom coefficients "%s".', customs)
        calibrator = cls.Calibrator(**spacecraft_coeffs)
        return calibrator

    @staticmethod
    def parse(value):
        """Cast lists to numpy arrays."""
        if isinstance(value, list):
            value = np.asarray(value)
        return value

    @classmethod
    def load_defaults(cls):
        """Read the default coefficients from file.

        Note:
            The pygac internal defaults are stored in data/calibration.json.
            The user can provide different defaults via the config file
            in the section "calibration" as option "coeffs_file", e.g.
            [calibration]
            coeffs_file = /path/to/user/default/coeffs.json
        """
        # check if the user has set a coefficient file
        config = pygac.configuration.get_config()
        coeffs_file = config.get("calibration", "coeffs_file", fallback='')
        # if no coeffs file has been specified, use the pygac defaults
        if not coeffs_file:
            coeffs_file = resource_filename('pygac', 'data/calibration.json')
        LOG.info('Use default coefficients from "%s"', coeffs_file)
        with open(coeffs_file, mode='r') as json_file:
            default_coeffs = json.load(json_file)
        # parse values in defaults
        cls.default_coeffs = {
            sat: {key: cls.parse(value) for key, value in defaults.items()}
            for sat, defaults in default_coeffs.items()
        }


def calibrate_solar(counts, chan, year, jday, spacecraft, corr=1, custom_coeffs=None):
    """Do the solar calibration and return reflectance (between 0 and 100).

    Arguments:
        counts (array) - counts for the given channel
        chan (array) - pygac internal channel index array
        year (int) - year
        jday (int) - day of year
        spacecraft (str) - pygac internal spacecraft name

    Optionals:
        corr (float) - reflectance correction multiplier (default = 1)
        custom_coeffs (dict) - custom calibration coefficients (default = None)

    Note:
        This function follows the solar calibration from PATMOS-x as described in
        Heidinger, A.K., W.C. Straka III, C.C. Molling, J.T. Sullivan, and X. Wu, (2010).
        Deriving an inter-sensor consistent calibration for the AVHRR solar reflectance data record.
        International Journal of Remote Sensing, 31:6493 - 6517.
    """
    # get the calibration coefficients for this spacecraft
    cal = Calibrator(spacecraft, custom_coeffs=custom_coeffs)

    # calculate the calibration slope with equation (6) in Heidinger et al 2010
    # S(t) = S_0*(100 + S_1*t + S_2*t^2) / 100,
    # where t is the time since launch expressed as years, S(t) is the calibration slope
    # and S_0, S_1 and S_2 are the coefficients of the quadratic fit.
    # See Fitting of Calibration Slope Equations in patmosx documentation
    # (CDRP-ATBD-0184 Rev. 2 03/29/2018 page 19)

    # Note that the this implementation does not take leap years into account!
    t = (year + jday / 365.0) - cal.l_date

    # apply slope equation for low and high gain coefficients
    stl = (cal.al[chan] * (100.0 + cal.bl[chan] * t
                           + cal.cl[chan] * t * t)) / 100.0
    sth = (cal.ah[chan] * (100.0 + cal.bh[chan] * t
                           + cal.ch[chan] * t * t)) / 100.0

    # Calculate the reflactance using equation (1) in Heidinger et al 2010
    # R = S*(C-D),
    # where R_cal is the value generated from the calibration and is referred to as a
    # scaled radiance, S is the calibration slope, C the measured count and D the dark count.
    # This equation is only valid for single gain instruments. Starting with the AVHRR/3 series
    # (1998, starting with NOAA-15, aka KLM), the channel-1, channel-2 and channel-3a require a
    # dual-gain calibration. See Appendix A in Heidinger et al 2010, for more information.
    # Long story short: The reflectance as function of instrument counts is a continuous piecewise
    # linear function of two line segments.
    #        R
    #        ^           *
    #        |          *             The reflactance R is given by
    #        |         *              R = S_l*(C - D)                 if C <= C_s
    # R(C_s) |--------*               R = R_cal(C_s) + S_h*(C - C_s)  if C > C_s
    #        |     *  |               where S_l and S_h are the low and high gain slopes,
    #        +--*-------------> C     D is the dark count, C_s is the gain switch.
    #     -D *       C_s
    # Note, that the implementation allows for an additional correction factor corr which defaults to one.
    if cal.c_s is not None:
        refl = np.where(counts <= cal.c_s[chan],
                        (counts - cal.c_dark[chan]) * stl * corr,
                        ((cal.c_s[chan] - cal.c_dark[chan]) * stl
                         + (counts - cal.c_s[chan]) * sth) * corr)
    else:
        refl = (counts - cal.c_dark[chan]) * stl * corr

    # Mask negative reflectances
    refl[refl < 0] = np.nan

    return refl


def calibrate_thermal(counts, prt, ict, space, line_numbers, channel, spacecraft, custom_coeffs=None):
    """Do the thermal calibration and return brightness temperatures (K).

    Arguments:
        counts (array) - counts for the given channel
        prt (array) - counts of the Platinum Resistance Thermometers (PRT)
        ict (array) - counts of the In-orbit Calibration Targets (ICT)
        space (array) - counts of cold space
        line_numbers (array) - line number index
        channel (array) - pygac internal channel index array
        spacecraft (str) - pygac internal spacecraft name

    Optionals:
        custom_coeffs (dict) - custom calibration coefficients (default = None)

    Note:
        This function follows steps 1 to 4 from the KLM guide section 7.1.2.4
        "Steps to Calibrate the AVHRR Thermal Channels"
    """
    # get the calibration coefficients for this spacecraft
    cal = Calibrator(spacecraft, custom_coeffs=custom_coeffs)

    # Shift channel index by three because from the possible channels [1, 2, 3a, 3b, 4, 5],
    # channels [3b, 4, 5] are the thermal channels
    chan = channel - 3

    lines, columns = counts.shape[:2]

    # Step 1. The temperature of the internal blackbody target is measured by four PRTs. In each
    # scanline, data words 18, 19 and 20 in the HRPT minor frame format contain three readings from
    # one of the four PRTs. (See Section 4.1.3) A different PRT is sampled each scanline; every fifth
    # scanline all three PRT values are set equal to 0 to indicate that a set of four PRTs has just been
    # sampled. The count value CPRT of each PRT is converted to temperature TPRT by the formula
    # T_PRT = d0 + d1*C_PRT + d2*C_PRT^2 + d3*C_PRT^3 + d4*C_PRT^4    (7.1.2.4-1)
    # The coefficients d0, d1, d2, d3 and d4 vary slightly for each PRT. Values for the coefficients are
    # found in Appendix D, in Table D.1-8 for NOAA-15 (coefficients d3 and d4 are 0 for NOAA-15),
    # Table D.2-9 for NOAA-16, Table D.3-3 for NOAA-17 and Table D.4-3 for NOAA-18. To
    # calculate the internal blackbody temperature TBB, NESDIS uses the simple average
    # T_BB = (T_PRT1 + T_PRT2 + T_PRT3 + T_PRT4)/4    (7.1.2.4-2)

    # Find the corresponding PRT values for a given line number
    # Why do we calculate this offset? Where does the threshold of prt_val < 50 come from?
    offset = 0

    for i, prt_val in enumerate(prt):
        if prt_val < 50:
            offset = i
            break

    # get the PRT index and fix some values (again based on threshold 50)
    iprt = (line_numbers - line_numbers[0] + 5 - offset) % 5

    ifix = np.where(np.logical_and(iprt == 1, prt < 50))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 1, prt > 50))
        prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])

    ifix = np.where(np.logical_and(iprt == 2, prt < 50))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 2, prt > 50))
        prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])

    ifix = np.where(np.logical_and(iprt == 3, prt < 50))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 3, prt > 50))
        prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])

    ifix = np.where(np.logical_and(iprt == 4, prt < 50))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 4, prt > 50))
        prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])

    # calculate temperature using equation (7.1.2.4-1)
    tprt = (cal.d[iprt, 0] + prt *
            (cal.d[iprt, 1] + prt *
             (cal.d[iprt, 2] + prt *
              (cal.d[iprt, 3] + prt *
               (cal.d[iprt, 4])))))

    # Note: the KLM Guide proposes to calculate the mean temperature using
    # equation (7.1.2.4-2). PyGAC follows a different Averaging approach.
    # It fills the zeros that mark a complete set of thermometer measurements
    # by interpolation, then it uses a weighting function (so far only equal
    # weighting) to convolve the temperatures (build a global average).
    # The same averaging technique is applied for ICTs and Space counts.
    zeros = iprt == 0
    nonzeros = np.logical_not(zeros)

    tprt[zeros] = np.interp((zeros).nonzero()[0],
                            (nonzeros).nonzero()[0],
                            tprt[nonzeros])

    if channel == 3:
        zeros = ict < 100
        nonzeros = np.logical_not(zeros)

        ict[zeros] = np.interp((zeros).nonzero()[0],
                               (nonzeros).nonzero()[0],
                               ict[nonzeros])
        zeros = space < 100
        nonzeros = np.logical_not(zeros)

        space[zeros] = np.interp((zeros).nonzero()[0],
                                 (nonzeros).nonzero()[0],
                                 space[nonzeros])

    # convolving and smoothing PRT, ICT and SPACE values
    if lines > 51:
        wlength = 51
    else:
        wlength = 3

    weighting_function = np.ones(wlength, dtype=float) / wlength
    tprt_convolved = np.convolve(tprt, weighting_function, 'same')
    ict_convolved = np.convolve(ict, weighting_function, 'same')
    space_convolved = np.convolve(space, weighting_function, 'same')

    # take care of the beginning and end
    tprt_convolved[0:(wlength - 1) // 2] = tprt_convolved[(wlength - 1) // 2]
    ict_convolved[0:(wlength - 1) // 2] = ict_convolved[(wlength - 1) // 2]
    space_convolved[0:(wlength - 1) // 2] = space_convolved[(wlength - 1) // 2]
    tprt_convolved[-(wlength - 1) // 2:] = tprt_convolved[-((wlength + 1) // 2)]
    ict_convolved[-(wlength - 1) // 2:] = ict_convolved[-((wlength + 1) // 2)]
    space_convolved[-(wlength - 1) // 2:] = \
        space_convolved[-((wlength + 1) // 2)]

    new_tprt = np.transpose(np.tile(tprt_convolved, (columns, 1)))
    new_ict = np.transpose(np.tile(ict_convolved, (columns, 1)))
    new_space = np.transpose(np.tile(space_convolved, (columns, 1)))

    # Step 2. The radiance NBB sensed in each thermal AVHRR channel from the internal blackbody
    # at temperature TBB is the weighted mean of the Planck function over the spectral response of the
    # channel. [...]. Each thermal channel has one equation, which uses a centroid wavenumber labmda_c and an
    # "effective" blackbody temperature TBB*. The two steps are:
    # TsBB = A + B*TBB    (7.1.2.4-3)
    # NBB = c1*nu_e^3/(exp(c2*nu_e/TsBB) - 1)    (7.1.2.4-3)
    # where c1 = 1.1910427e-5 mW/m^2/sr/cm^{-4}, c2 = 1.4387752 cm K

    # calibrating thermal channel

    tBB = new_tprt
    tsBB = cal.a[chan] + cal.b[chan] * tBB
    nBB_num = (1.1910427 * 0.000010) * cal.c_wn[chan] ** 3
    nBB = nBB_num / (np.exp((1.4387752 * cal.c_wn[chan]) / tsBB) - 1.0)

    # Step 3. Output from the two in-orbit calibration targets is used to compute a linear estimate of
    # the Earth scene radiance NE. Each scanline, the AVHRR views the internal blackbody target and
    # outputs 10 count values for each of the three thermal channel detectors; these are found in words
    # 23 to 52 in the HRPT data stream. When the AVHRR views cold space, 10 counts from each of
    # the five channel sensors are output and placed into words 52 to 102. (Table 4.1.3-1 describes
    # how these data are multiplexed.) Count values for each channel are averaged together to smooth
    # our random noise; often counts from five consecutive scanlines are averaged because it takes five
    # lines to obtain a set of all four PRT measurements. The average blackbody count CBB and the
    # average space count CS, together with blackbody radiance NBB and space radiance NS, explained
    # in the next paragraph, are used to compute the linear radiance estimate NLIN,
    # NLIN = NS + (NBB - NS)*(CS - CE)/(CS -CBB)    (7.1.2.4-5)
    # where CE is the AVHRR count output when it views one of the 2,048 Earth targets.
    # The Mercury-Cadmium-Telluride detectors used for channels 4 and 5 have a nonlinear response
    # to incoming radiance. Pre-launch laboratory measurements show that:
    #     a. scene radiance is a slightly nonlinear (quadratic) function of AVHRR output count,
    #     b. the nonlinearity depends on the AVHRR operating temperature.
    # It is assumed that the nonlinear response will persist in orbit. For the NOAA KLM series of
    # satellites, NESDIS uses a radiance-based nonlinear correction method. In this method, the linear
    # radiance estimate is first computed using a non-zero radiance of space, the NS term in Equation
    # 7.1.2.4-5. Then, the linear radiance value is input into a quadratic equation to generate the
    # nonlinear radiance correction NCOR:
    # NCOR = b0 + b1*NLIN + b2*NLIN^2    (7.1.2.4-6)
    # Finally, the Earth scene radiance is obtained by adding NCOR to NLIN
    # NE = NLIN + NCOR

    # Note: PyGAC increases the linear correction coefficient by 1 to accout for the linear part!
    # Therefore, NE = NCOR in PyGAC

    Nlin = (cal.n_s[chan] +
            (((nBB - cal.n_s[chan])
              * (new_space - counts.astype(float)))
             / (new_space - new_ict)))
    Ncor = cal.b0[chan] + Nlin * (cal.b1[chan] + cal.b2[chan] * Nlin)
    Ne = Ncor

    # Step 4. Data users often convert the computed Earth scene radiance value NE into an equivalent
    # blackbody temperature TE. This temperature is defined by simple inverting the steps used to
    # calculate the radiance NE sensed by an AVHRR channel from an emitting blackbody at a
    # temperature TE. The two-step process is:
    # TsE = c2*nu_c / ln(1 + (c1*nu_c/NE))    (7.1.2.4-8)
    # TE = (TsE - A)/B    (7.1.2.4-9)

    tsE = ((1.4387752 * cal.c_wn[chan])
           / np.log(1.0 + nBB_num / Ne))
    bt = (tsE - cal.a[chan]) / cal.b[chan]

    # Why do we do this on channel 3b?
    if chan == 0:
        bt = np.where((counts - new_space) >= 0, 0.0, bt)

    # Mask values outside valid range
    bt = np.where(np.logical_or(bt < 170.0, bt > 350.0), np.nan, bt)

    return bt
