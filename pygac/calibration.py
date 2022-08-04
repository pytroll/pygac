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
from enum import Enum
import sys
import logging
import numpy as np
import json
import hashlib
import warnings
import datetime as dt
from collections import namedtuple
from pkg_resources import resource_filename


LOG = logging.getLogger(__name__)


class CoeffStatus(Enum):
    """Indicates the status of calibration coefficients."""
    NOMINAL = 'nominal'
    PROVISIONAL = 'provisional'
    EXPERIMENTAL = 'experimental'


class Calibrator(object):
    """Factory class to create namedtuples holding the calibration coefficients.

    Attributes:
        fields: coefficient names
        Calibrator: namedtuple constructor
        default_coeffs: dictonary containing default values for all spacecrafts
    """
    version_hashs = {
        '963af9b66268475ed500ad7b37da33c5': {
            'name': 'PATMOS-x, v2017r1',
            'status': CoeffStatus.NOMINAL
        },
        '689386c822de18a07194ac7fd71652ea': {
            'name': 'PATMOS-x, v2017r1, with provisional coefficients for MetOp-C',
            'status': CoeffStatus.PROVISIONAL
        }
    }
    fields = [
        "dark_count", "gain_switch", "s0", "s1", "s2", "b",  # "b0", "b1", "b2",
        "centroid_wavenumber", "space_radiance", "to_eff_blackbody_intercept",
        "to_eff_blackbody_slope", "date_of_launch", "d", "spacecraft", "version"
    ]

    Calibrator = namedtuple('Calibrator', fields)
    default_coeffs = None
    default_file = None
    default_version = None

    def __new__(cls, spacecraft, custom_coeffs=None, coeffs_file=None):
        """Creates a namedtuple for calibration coefficients of a given spacecraft

        Args:
            spacecraft (str): spacecraft name in pygac convention
            custom_coeffs (dict): custom coefficients (optional)
            coeffs_file (str): path to coefficents file (optional)

        Returns:
            calibrator (namedtuple): calibration coefficients

        Note:
            The coefficients in coeffs_file serve as default values if no custom_coeffs
            are given. If omitted, Calibrator uses the PyGAC internal defaults.
        """
        if cls.default_coeffs is None or cls.default_file != coeffs_file:
            cls.default_file = coeffs_file
            cls.default_coeffs, cls.default_version = cls.read_coeffs(coeffs_file)
        if custom_coeffs:
            LOG.info('Using following custom coefficients "%s".', custom_coeffs)
        customs = custom_coeffs or {}
        defaults = cls.default_coeffs[spacecraft]
        coeffs = defaults.copy()
        coeffs.update(customs)

        # transpose the coefficient order from channel - coeff to coeff - channel
        # and store as arrays for vectorized calls of calibration functions
        arraycoeffs = dict.fromkeys(cls.fields)
        # visible channels
        for key in ("dark_count", "gain_switch", "s0", "s1", "s2"):
            arraycoeffs[key] = np.array([
                coeffs[channel][key]
                for channel in ('channel_1', 'channel_2', 'channel_3a')
            ], dtype=float)
        # thermal channels
        for key in ("centroid_wavenumber", "space_radiance",
                    "to_eff_blackbody_intercept", "to_eff_blackbody_slope"):
            arraycoeffs[key] = np.array([
                coeffs[channel][key]
                for channel in ('channel_3b', 'channel_4', 'channel_5')
            ], dtype=float)
        arraycoeffs["b"] = np.array([
            [
                coeffs[channel][key]
                for key in ("b0", "b1", "b2")
            ]
            for channel in ('channel_3b', 'channel_4', 'channel_5')
        ], dtype=float)
        # thermometers
        # Note, that "thermometer_0" does not exists, and is filled with zeros to
        # account for the PRT reset every fifth scanline
        arraycoeffs["d"] = np.array([
            [
                coeffs.get("thermometer_{0}".format(t), {}).get("d{0}".format(d), 0.0)
                for t in range(5)
            ]
            for d in range(5)
        ], dtype=float)
        # parse date of launch
        date_of_launch_str = coeffs["date_of_launch"].replace('Z', '+00:00')
        if sys.version_info < (3, 7):
            # Note that here any time information is lost
            import dateutil.parser
            date_of_launch = dateutil.parser.parse(date_of_launch_str)
        else:
            # datetime.fromisoformat() was introduced in Python-3.7
            date_of_launch = dt.datetime.fromisoformat(
                date_of_launch_str).astimezone(dt.timezone.utc)
        # remove time zone information (easier to handle in calculations)
        arraycoeffs["date_of_launch"] = date_of_launch.replace(tzinfo=None)
        arraycoeffs["spacecraft"] = spacecraft
        arraycoeffs["version"] = cls.default_version
        if custom_coeffs:
            arraycoeffs["version"] = None
        # create namedtuple
        calibrator = cls.Calibrator(**arraycoeffs)
        return calibrator

    @classmethod
    def read_coeffs(cls, coeffs_file):
        """Read calibration coefficients for all satellites from file.

        Args:
            coeffs_file (str): path to coefficients file

        Returns:
            coeffs (dict): dictionary containing coefficients for all satellites
            version (str): version of the coefficients (None if unknown)
        """
        if coeffs_file:
            LOG.info('Read calibration coefficients from "%s"', coeffs_file)
        else:
            LOG.debug("Read PyGAC internal calibration coefficients.")
            coeffs_file = resource_filename('pygac', 'data/calibration.json')
        with open(coeffs_file, mode='rb') as json_file:
            content = json_file.read()
            coeffs = json.loads(content)
            version = cls._get_coeffs_version(content)
        return coeffs, version

    @classmethod
    def _get_coeffs_version(cls, coeff_file_content):
        """Determine coefficient version."""
        md5_hash = hashlib.md5(coeff_file_content)
        digest = md5_hash.hexdigest()
        version_dict = cls.version_hashs.get(
            digest,
            {'name': None, 'status': None}
        )
        version = version_dict['name']
        status = version_dict['status']
        if version is None:
            warning = "Unknown calibration coefficients version!"
            warnings.warn(warning, RuntimeWarning)
            LOG.warning(warning)
        else:
            LOG.info('Identified calibration coefficients version "%s".',
                     version)
            if status != CoeffStatus.NOMINAL:
                warning = 'Using {} calibration coefficients'.format(status)
                warnings.warn(warning, RuntimeWarning)
                LOG.warning(warning)
        return version

    @staticmethod
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


def calibrate_solar(counts, chan, year, jday, cal, corr=1):
    """Do the solar calibration and return scaled radiance.

    Arguments:
        counts (array) - raw counts for the given channels (options 1, 2[, 3A if available & active])
        chan (array) - pygac internal channel index array
        year (int) - year
        jday (int) - day of year
        cal (namedtuple) - spacecraft specific calibration coefficients, see Calibrator

    Optional:
        corr (float) - depricated - reflectance correction multiplier (default = 1)

    Returns:
        r_cal (array) - scaled radiance

    Note:
        This function and documentation follows the time-dependent solar calibration as described in:
        Heidinger, A.K., W.C. Straka III, C.C. Molling, J.T. Sullivan, and X. Wu, (2010).
        "Deriving an inter-sensor consistent calibration for the AVHRR solar reflectance data record.",
        International Journal of Remote Sensing, 31:6493 - 6517.
    """
    # Step 1. Obtain the calibration slope using equation (6) in Heidinger et al 2010
    # S(t) = S_0*(100 + S_1*t + S_2*t^2) / 100,
    # where t is the time since launch expressed as years, S(t) is the calibration slope
    # and S_0, S_1 and S_2 are the coefficients of the quadratic fit.
    # See also section "Fitting of Calibration Slope Equations" in PATMOS-x documentation
    # (CDRP-ATBD-0184 Rev. 2 03/29/2018 page 19)

    # Note that the this implementation does not take leap years into account!
    # Using datetime objects would include it automatically
    # sensing_date = dt.strptime('{0}.{1}'.format(year, jday), '%Y.%j')
    # delta = sensing_date - cal.date_of_launch
    # t = delta.total_seconds() / 31557600  # divided by seconds of a Julian year
    l_date = Calibrator.date2float(cal.date_of_launch)
    t = (year + jday / 365.0) - l_date

    # Note: splitting the calibration slope is needed to reproduce old results and may disappear in future,
    #       because actually there is only one set of slope parameters defined for single-gain counts
    #       as described in Heidinger et al. 2010. See Step 2 for more information.
    # Note that in case of a single-gain instrument, all gain_switch parameters are set to NaN.
    if np.isnan(cal.gain_switch).all():
        glow = ghigh = np.ones(3)
    else:
        glow = np.array([0.5, 0.5, 0.25])
        ghigh = np.array([1.5, 1.5, 1.75])
    # especially the rounding to three digits is crutial to exectly reproduce the original PATMOS-x values.
    al, bl, cl = np.round(glow*cal.s0, 3), cal.s1, cal.s2
    ah, bh, ch = np.round(ghigh*cal.s0, 3), cal.s1, cal.s2

    # apply slope equation for low and high gain coefficients
    stl = (al[chan] * (100.0 + bl[chan] * t + cl[chan] * t * t)) / 100.0
    sth = (ah[chan] * (100.0 + bh[chan] * t + ch[chan] * t * t)) / 100.0

    # Step 2. Calculate the scaled radiance using equation (1) in Heidinger et al 2010
    # R_cal = S*(C-D),
    # where R_cal is the value generated from the calibration and is referred to as a
    # scaled radiance, S is the calibration slope, C the measured count and D the dark count.
    # This equation is only valid for single-gain instruments. Starting with the AVHRR/3 series
    # (from NOAA-15 onwards), the channel-1, 2 and 3a require a dual-gain calibration.
    # The conversion for channel-1 and 2 is given in the appendix, equation (A1) and (A2).
    # In general, these equations can be written as
    # C(C_dg) = D + G_low*(C_dg-D),            if C_dg <= B_dg
    # C(C_dg) = C(B_dg) + G_high*(C_dg-B_dg),  otherwise
    # where C_dg is the measured dual-gain counts, B_dg is the dual-gain switch and G_low/high
    # are the gain factors for the low and high count region.
    # Quote from the book "Remote Sensing Time Series: Revealing Land Surface Dynamics":
    # > Another change in design of the AVHRR/3 instrument was the introduction of a dual-gain feature
    # > for the reflective channels 1, 2 and 3A. In order to improve the radiometric resolution of the
    # > instrument for low reflectance targets, the dynamic range of the instrument was divided equally
    # > in two ranges, i.e. nominally from 0 to 500 counts and from 500 to 1,000 counts. For channels 1
    # > and 2 half of the available Digital Number (DN) range is assigned to the low albedo range from
    # > 0 to 25% with the other half to the high albedo range from 26 to 100%. This allows for an increase
    # > in the radiometric resolution for dark targets. For channel 3A, the split between low and high
    # > albedo range is set at 12.5% albedo (Rao and Sullivan 2001).
    # The gain factors are given by the ratio of the fraction of albedo range to the fraction of count range
    # for the given count region. From the information given by the book quote, we get
    # G_low = 25% / 50% = 0.5 for channel-1 and 2
    # G_low = 12.5% / 50% = 0.25 for channel-3a
    # G_high = (100% - 25%) / (100% - 50%) = 1.5 for channel-1 and 2
    # G_high = (100% - 12.5%) / (100% - 50%) = 1.75 for channel-3a
    # Inserting the converted dual-gain counts into equation (1) yields the scaled radiance equation, which
    # is a continuous piecewise linear function of two line segments.
    #            R_cal
    #            ^           *         R_cal = S*G_low*(C_dg-D),                        if C_dg <= B_dg
    #            |          *          R_cal = S*G_low*(B_dg-D) + S*G_high*(C_dg-B_dg), otherwise
    #            |         *
    # R_cal(B_dg)|--------*
    #            |     *  |
    #          0 +--*-------------> C_dg
    #            *  D    B_dg
    # Note, that in the former implementation, there was a distinction beteen low and high gain slopes
    # given by S_low/high = S*G_low/high. which only affects S0 in equation (6) in Heidinger et al 2010.
    # Furthermore, the implementation allows for an additional correction factor corr which defaults to 1. (depricated)
    d = cal.dark_count[chan]
    b_dg = cal.gain_switch[chan]
    # Note that in case of a single-gain instrument, all gain_switch parameters are set to NaN.
    if not np.isnan(cal.gain_switch).all():
        r_cal = np.where(
            counts <= b_dg,
            (counts - d)*stl,
            (b_dg - d)*stl + (counts - b_dg)*sth
        )
    else:
        r_cal = stl*(counts - d)
    # apply depricated correction
    if corr != 1:
        warnings.warn(
            "Using the 'corr' argument is depricated in favor of making the units"
            " of the function result clear. Please make any unit conversion outside this function.",
            DeprecationWarning
        )
        r_cal *= corr

    # Mask negative scaled radiances
    r_cal[r_cal < 0] = np.nan

    return r_cal


def calibrate_thermal(counts, prt, ict, space, line_numbers, channel, cal):
    """Do the thermal calibration and return brightness temperatures (K).

    Arguments:
        counts (array) - counts for the given channel (options: 3B (if active), 4, 5)
        prt (array) - counts of the Platinum Resistance Thermometers (PRT)
        ict (array) - counts of the In-orbit Calibration Targets (ICT)
        space (array) - counts of cold space
        line_numbers (array) - line number index
        channel (array) - pygac internal channel index array
        cal (namedtuple) - spacecraft specific calibration coefficients, see Calibrator

    Note:
        This function and documentation follows steps 1 to 4 from the  KLM User's Guide
        (Robel, J. (2009). NOAA KLM user's guide with NOAA-N,-P supplement. NOAA KLM Users
        Guide - August 2014 Revision) section 7.1.2.4 "Steps to Calibrate the AVHRR Thermal Channels",
        and the smoothing approach by Trishchenko (2002).
        The correction method for the non-linear response of the Mercury-Cadmium-Telluride detectors used
        for channels 4 and 5 is based on Walton et al. (1998)
    """
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
    # Tprt = d0 + d1*Cprt + d2*Cprt^2 + d3*Cprt^3 + d4*Cprt^4
    # Note: First dimension of cal.d are the five coefficient indicees
    tprt = np.polynomial.polynomial.polyval(prt, cal.d[:, iprt], tensor=False)

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

    # Thresholds to flag missing/wrong data for interpolation
    ict_threshold = 100
    space_threshold = 100
    if channel == 3:
        zeros = ict < ict_threshold
        nonzeros = np.logical_not(zeros)

        ict[zeros] = np.interp((zeros).nonzero()[0],
                               (nonzeros).nonzero()[0],
                               ict[nonzeros])
        zeros = space < space_threshold
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
    # where the constants of the Planck function are defined as c1 = 2*h*c^2, c2 = h*c/k_B.
    # constatns
    c1 = 1.1910427e-5  # mW/m^2/sr/cm^{-4}
    c2 = 1.4387752  # cm K
    # coefficients
    A = cal.to_eff_blackbody_intercept[chan]
    B = cal.to_eff_blackbody_slope[chan]
    nu_c = cal.centroid_wavenumber[chan]
    nS = cal.space_radiance[chan]
    b = cal.b[chan, 0:3]  # the second index are the three polynomial coefficients
    # variables
    tBB = new_tprt
    cS = new_space
    cE = counts.astype(float)
    cBB = new_ict

    tsBB = A + B*tBB
    nBB_num = c1 * nu_c**3
    nBB = nBB_num / (np.exp((c2 * nu_c) / tsBB) - 1.0)

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
    # NLIN = NS + (NBB - NS)*(CS - CE)/(CS - CBB)    (7.1.2.4-5)
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

    # Note: For channel 3B, the non-linear correction coefficients are set to zero to use the same equation.
    Nlin = nS + (nBB - nS)*(cS - cE)/(cS - cBB)
    Ncor = np.polynomial.polynomial.polyval(Nlin, b[0:3], tensor=False)
    Ne = Nlin + Ncor

    # Step 4. Data users often convert the computed Earth scene radiance value NE into an equivalent
    # blackbody temperature TE. This temperature is defined by simple inverting the steps used to
    # calculate the radiance NE sensed by an AVHRR channel from an emitting blackbody at a
    # temperature TE. The two-step process is:
    # TsE = c2*nu_c / ln(1 + (c1*nu_c^3/NE))    (7.1.2.4-8)
    # TE = (TsE - A)/B    (7.1.2.4-9)
    tsE = c2*nu_c / np.log(1.0 + nBB_num / Ne)
    bt = (tsE - A) / B

    # Why do we do this on channel 3b?
    if chan == 0:
        bt = np.where((counts - new_space) >= 0, 0.0, bt)

    # Mask values outside valid range
    bt = np.where(np.logical_or(bt < 170.0, bt > 350.0), np.nan, bt)

    return bt
