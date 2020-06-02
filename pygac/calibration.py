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
import warnings
import datetime as dt
import dateutil.parser
from collections import namedtuple, defaultdict
from pkg_resources import resource_filename


import pygac.configuration

LOG = logging.getLogger(__name__)


def date2float(date, decimals=5):
    """Convert date to year float.

    Argument
        date (datetime.datetime) - date
        decimals (int or None) - rounding precision if None, do not round (default=5)

    Return
        date_float (float) - date as year

    Note
        rounding to the 5th decimal reproduces the original float values from patmos-x

    Example:
        date2float('2000-07-02') == 2000.5
        because the 2nd of July was the middle day of the leap year 2000
    """
    year = date.year
    days_in_year = (dt.datetime(year+1, 1, 1) - dt.datetime(year, 1, 1)).days
    diff = date - dt.datetime(year, 1, 1)
    seconds = diff.total_seconds()
    date_float = date.year + seconds/(days_in_year*24*3600)
    if decimals is not None:
        date_float = round(date_float, decimals)
    return date_float


def float2date(date_float):
    """Convert date float to date.

    Argument
        date_float (float) - date as year

    Return
        date (datetime.datetime) - date

    Note
        This is the reverse function of date2float.
    """
    year = int(date_float)
    days_in_year = (dt.datetime(year+1, 1, 1) - dt.datetime(year, 1, 1)).days
    seconds = date_float*days_in_year*24*3600 - year*days_in_year*24*3600
    diff = dt.timedelta(seconds=seconds)
    date = dt.datetime(year, 1, 1) + diff
    return date


def new2old_coeffs(new_coeffs):
    """convert new coefficients to old coefficients"""
    old_coeffs = defaultdict(list)
    date_of_launch = dateutil.parser.parse(new_coeffs['date_of_launch'])
    old_coeffs['l_date'] = date2float(date_of_launch)
    for i, ch in enumerate(['1', '2', '3a']):
        for slope, char in enumerate('abc'):
            for gain in ['high', 'low']:
                old_coeffs['{0}{1}'.format(char, gain[0])].append(
                    new_coeffs['channel_{0}'.format(ch)]['gain_{0}_s{1}'.format(gain, slope)]
                )
        old_coeffs['c_dark'].append(new_coeffs['channel_{0}'.format(ch)]['dark_count'])
        if new_coeffs['channel_1'].get('gain_switch') is not None:
            old_coeffs['c_s'].append(new_coeffs['channel_{0}'.format(ch)]['gain_switch'])
    old_coeffs['d'].append(5*[0.0])
    for prt in range(1, 5):
        old_coeffs['d'].append([
            new_coeffs['thermometer_{0}'.format(prt)]['c{0}'.format(i)]
            for i in range(5)
        ])
    for ch in ['3b', '4', '5']:
        old_coeffs['n_s'].append(new_coeffs['channel_{0}'.format(ch)]['space_radiance'])
        old_coeffs['c_wn'].append(new_coeffs['channel_{0}'.format(ch)]['centroid_wavenumber'])
        old_coeffs['a'].append(new_coeffs['channel_{0}'.format(ch)]['to_eff_blackbody_intercept'])
        old_coeffs['b'].append(new_coeffs['channel_{0}'.format(ch)]['to_eff_blackbody_slope'])
        for j in range(3):
            # revert to original definition of the linear coefficient (not including +1)
            # if j == 1:
            #    old_coeffs['b{0}'.format(j)].append(1 + new_coeffs['channel_{0}'.format(ch)][
            #        'radiance_correction_c{0}'.format(j)])
            # else:
            #    old_coeffs['b{0}'.format(j)].append(new_coeffs['channel_{0}'.format(ch)][
            #        'radiance_correction_c{0}'.format(j)])
            old_coeffs['b{0}'.format(j)].append(new_coeffs['channel_{0}'.format(ch)][
                'radiance_correction_c{0}'.format(j)])
    return dict(old_coeffs)


def old2new_coeffs(old_coeffs):
    """convert old coefficients to new coefficients"""
    new_coeffs = {}
    for i, ch in enumerate(['1', '2', '3a']):
        new_coeffs['channel_{0}'.format(ch)] = {}
        for slope, char in enumerate('abc'):
            for gain in ['high', 'low']:
                new_coeffs['channel_{0}'.format(ch)]['gain_{0}_s{1}'.format(gain, slope)] = old_coeffs[
                    '{0}{1}'.format(char, gain[0])][i]
        new_coeffs['channel_{0}'.format(ch)]['dark_count'] = old_coeffs['c_dark'][i]
        new_coeffs['channel_{0}'.format(ch)]['gain_switch'] = old_coeffs.get('c_s', 3*[None])[i]
    new_coeffs['date_of_launch'] = str(float2date(old_coeffs['l_date']))
    for prt in range(1, 5):
        new_coeffs['thermometer_{0}'.format(prt)] = {}
        for i in range(5):
            new_coeffs['thermometer_{0}'.format(prt)]['c{0}'.format(i)] = old_coeffs['d'][prt][i]
    for i, ch in enumerate(['3b', '4', '5']):
        new_coeffs['channel_{0}'.format(ch)] = {}
        new_coeffs['channel_{0}'.format(ch)]['space_radiance'] = old_coeffs['n_s'][i]
        new_coeffs['channel_{0}'.format(ch)]['centroid_wavenumber'] = old_coeffs['c_wn'][i]
        new_coeffs['channel_{0}'.format(ch)]['to_eff_blackbody_intercept'] = old_coeffs['a'][i]
        new_coeffs['channel_{0}'.format(ch)]['to_eff_blackbody_slope'] = old_coeffs['b'][i]
        for j in range(3):
            new_coeffs['channel_{0}'.format(ch)]['radiance_correction_c{0}'.format(j)] = old_coeffs[
                'b{0}'.format(j)][i]
            # revert to original definition of the linear coefficient (not including +1)
            # if j == 1:
            #    new_coeffs['channel_{0}'.format(ch)]['radiance_correction_c{0}'.format(j)] -= 1
    return new_coeffs


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
        spacecraft_coeffs = {}  # dict.fromkeys(cls.fields)
        spacecraft_coeffs.update(defaults)
        spacecraft_coeffs.update(customs)
        if custom_coeffs:
            LOG.info('Using following custom coefficients "%s".', customs)
        # Note: The follwoing conversion from new to old coefficient names should
        #       disappear in a future version.
        #       The commit intorducing this comment contains everything to convert
        #       new to old names and vice versa
        _spacecraft_coeffs = dict.fromkeys(cls.fields)
        _spacecraft_coeffs.update({
            key: cls.parse(value)  # parse lists to arrays
            for key, value in new2old_coeffs(spacecraft_coeffs).items()
        })
        calibrator = cls.Calibrator(**_spacecraft_coeffs)
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
            cls.default_coeffs = json.load(json_file)


def calibrate_solar(counts, chan, year, jday, spacecraft, corr=1, custom_coeffs=None):
    """Do the solar calibration and return scaled radiance.

    Arguments:
        counts (array) - raw counts for the given channels (options 1, 2[, 3A if available & active])
        chan (array) - pygac internal channel index array
        year (int) - year
        jday (int) - day of year
        spacecraft (str) - pygac internal spacecraft name

    Optionals:
        corr (float) - depricated - reflectance correction multiplier (default = 1)
        custom_coeffs (dict) - custom calibration coefficients (default = None)

    Returns:
        r_cal (array) - scaled radiance

    Note:
        This function and documentation follows the time-dependent solar calibration as described in:
        Heidinger, A.K., W.C. Straka III, C.C. Molling, J.T. Sullivan, and X. Wu, (2010).
        "Deriving an inter-sensor consistent calibration for the AVHRR solar reflectance data record.",
        International Journal of Remote Sensing, 31:6493 - 6517.
    """
    if corr != 1:
        warnings.warn(
            "Using the 'corr' argument is depricated in favor of making the units"
            " of the function result clear. Please make any unit conversion outside this function.",
            DeprecationWarning
        )

    # get the calibration coefficients for this spacecraft
    cal = Calibrator(spacecraft, custom_coeffs=custom_coeffs)

    # Step 1. Obtain the calibration slope using equation (6) in Heidinger et al 2010
    # S(t) = S_0*(100 + S_1*t + S_2*t^2) / 100,
    # where t is the time since launch expressed as years, S(t) is the calibration slope
    # and S_0, S_1 and S_2 are the coefficients of the quadratic fit.
    # See section "Fitting of Calibration Slope Equations" in PATMOS-x documentation
    # (CDRP-ATBD-0184 Rev. 2 03/29/2018 page 19)

    # Note that the this implementation does not take leap years into account!
    t = (year + jday / 365.0) - cal.l_date

    # apply slope equation for low and high gain coefficients
    stl = (cal.al[chan] * (100.0 + cal.bl[chan] * t
                           + cal.cl[chan] * t * t)) / 100.0
    sth = (cal.ah[chan] * (100.0 + cal.bh[chan] * t
                           + cal.ch[chan] * t * t)) / 100.0

    # Step 2. Calculate the scaled radiance using equation (1) in Heidinger et al 2010
    # R_cal = S*(C-D),
    # where R_cal is the value generated from the calibration and is referred to as a
    # scaled radiance, S is the calibration slope, C the measured count and D the dark count.
    # This equation is only valid for single gain instruments. Starting with the AVHRR/3 series
    # (from NOAA-15 onwards), the channel-1, 2 and 3a require a dual-gain calibration. The dual-gain
    # counts have been converted to equivalent single-gain counts determined by a gain switch to be
    # linearly propotional to radiance (Appendix A in Heidinger et al, 2010).
    # Long story short: The scaled radiance as function of instrument counts is a continuous piecewise
    # linear function of two line segments.
    #            R_cal
    #            ^           *
    #            |          *             The scaled radiance R_cal is given by
    #            |         *              R_cal = S_l*(C - D)                 if C <= C_s
    # R_cal(C_s) |--------*               R_cal = R_cal(C_s) + S_h*(C - C_s)  if C > C_s
    #            |     *  |               where S_l and S_h are the low and high gain slopes,
    #          0 +--*-------------> C     D is the dark count, C_s is the gain switch.
    #            *  D    C_s
    # Note, that the implementation allows for an additional correction factor corr which defaults to one. (depricated)
    if cal.c_s is not None:
        r_cal = np.where(counts <= cal.c_s[chan],
                        (counts - cal.c_dark[chan]) * stl * corr,
                        ((cal.c_s[chan] - cal.c_dark[chan]) * stl
                         + (counts - cal.c_s[chan]) * sth) * corr)
    else:
        r_cal = (counts - cal.c_dark[chan]) * stl * corr

    # Mask negative scaled radiances
    r_cal[r_cal < 0] = np.nan

    return r_cal


def calibrate_thermal(counts, prt, ict, space, line_numbers, channel, spacecraft, custom_coeffs=None):
    """Do the thermal calibration and return brightness temperatures (K).

    Arguments:
        counts (array) - counts for the given channel (options: 3B (if active), 4, 5)
        prt (array) - counts of the Platinum Resistance Thermometers (PRT)
        ict (array) - counts of the In-orbit Calibration Targets (ICT)
        space (array) - counts of cold space
        line_numbers (array) - line number index
        channel (array) - pygac internal channel index array
        spacecraft (str) - pygac internal spacecraft name

    Optionals:
        custom_coeffs (dict) - custom calibration coefficients (default = None)

    Note:
        This function and documentation follows steps 1 to 4 from the  KLM User’s Guide 
        (Robel, J. (2009). NOAA KLM user's guide with NOAA-N,-P supplement. NOAA KLM Users 
        Guide –August 2014 Revision) section 7.1.2.4 "Steps to Calibrate the AVHRR Thermal Channels",
        and the smoothing approach by Trishchenko (2002).
        The correction method for the non-linear response of the Mercury-Cadmium-Telluride detectors used
        for channels 4 and 5 is based on Walton et al. (1998)
    """
    # get the calibration coefficients for this spacecraft
    cal = Calibrator(spacecraft, custom_coeffs=custom_coeffs)

    # Shift channel index by three to obtain thermal channels [3b, 4, 5].
    chan = channel - 3

    lines, columns = counts.shape[:2]

    # Step 1. The temperature of the internal blackbody target is measured by four platinum resistance
    # thermometers (PRT)s. In each scanline, data words 18, 19 and 20 in the HRPT minor frame format contain
    # three readings from one of the four PRTs. (See Section 4.1.3) A different PRT is sampled each scanline;
    # every fifth scanline all three PRT values are set equal to 0 to indicate that a set of four PRTs has
    # just been sampled. The count value CPRT of each PRT is converted to temperature TPRT by the formula
    # T_PRT = d0 + d1*C_PRT + d2*C_PRT^2 + d3*C_PRT^3 + d4*C_PRT^4    (7.1.2.4-1)
    # The coefficients d0, d1, d2, d3 and d4 vary slightly for each PRT. Values for the coefficients are
    # found in Appendix D, in Table D.1-8 for NOAA-15 (coefficients d3 and d4 are 0 for NOAA-15),
    # Table D.2-9 for NOAA-16, Table D.3-3 for NOAA-17 and Table D.4-3 for NOAA-18. To
    # calculate the internal blackbody temperature TBB, NESDIS uses the simple average
    # T_BB = (T_PRT1 + T_PRT2 + T_PRT3 + T_PRT4)/4    (7.1.2.4-2)

    # Find the corresponding PRT values for a given line number
    # Note that the prt values are the average value of the three readings from one of the four
    # PRTs. See reader.get_telemetry implementations.
    prt_threshold = 50  # empirically found and set by Abhay Devasthale
    offset = 0

    for i, prt_val in enumerate(prt):
        # According to the KLM Guide the fill value between PRT measurments is 0, but we search
        # for the first measurment gap using the threshold. Is this on purpose?
        if prt_val < prt_threshold:
            offset = i
            break

    # get the PRT index, iprt equals to 0 corresponds to the measurement gaps
    iprt = (line_numbers - line_numbers[0] + 5 - offset) % 5

    # fill measured values below threshold by interpolation
    ifix = np.where(np.logical_and(iprt == 1, prt < prt_threshold))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 1, prt > prt_threshold))
        prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])

    ifix = np.where(np.logical_and(iprt == 2, prt < prt_threshold))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 2, prt > prt_threshold))
        prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])

    ifix = np.where(np.logical_and(iprt == 3, prt < prt_threshold))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 3, prt > prt_threshold))
        prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])

    ifix = np.where(np.logical_and(iprt == 4, prt < prt_threshold))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 4, prt > prt_threshold))
        prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])

    # calculate PRT temperature using equation (7.1.2.4-1) KLM Guide
    tprt = (cal.d[iprt, 0] + prt *
            (cal.d[iprt, 1] + prt *
             (cal.d[iprt, 2] + prt *
              (cal.d[iprt, 3] + prt *
               (cal.d[iprt, 4])))))

    # Note: the KLM Guide proposes to calculate the mean temperature using equation (7.1.2.4-2).
    # PyGAC follows the smoothing approach by Trishchenko (2002), i.e.
    # filling the zeros that mark a complete set of thermometer measurements
    # by interpolation, and then using a weighting function (so far only equal
    # weighting) to convolve the temperatures to calculate a moving average of a given window size.
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
        wlength = 51  # empirically found and set by Abhay Devasthale
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
    space_convolved[-(wlength - 1) // 2:] = space_convolved[-((wlength + 1) // 2)]

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
    # where CE is the AVHRR count output when it views one of the reference Earth targets.
    # While the detector in channel 3B has a linear response, the Mercury-Cadmium-Telluride detectors used
    # for channels 4 and 5 have a nonlinear response to incoming radiance. Pre-launch laboratory measurements show that:
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

    # Note: For channel 3B, the correction coefficients are set to zero to use the same equation.
    Nlin = (cal.n_s[chan] +
            (((nBB - cal.n_s[chan])
              * (new_space - counts.astype(float)))
             / (new_space - new_ict)))
    Ncor = cal.b0[chan] + Nlin * (cal.b1[chan] + cal.b2[chan] * Nlin)
    Ne = Nlin + Ncor

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
