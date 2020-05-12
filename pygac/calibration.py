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
    """Do the solar calibration and return reflectance (between 0 and 100)."""
    cal = Calibrator(spacecraft, custom_coeffs=custom_coeffs)

    t = (year + jday / 365.0) - cal.l_date
    stl = (cal.al[chan] * (100.0 + cal.bl[chan] * t
                           + cal.cl[chan] * t * t)) / 100.0
    sth = (cal.ah[chan] * (100.0 + cal.bh[chan] * t
                           + cal.ch[chan] * t * t)) / 100.0
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
    """Do the thermal calibration and return brightness temperatures (K)."""
    cal = Calibrator(spacecraft, custom_coeffs=custom_coeffs)

    chan = channel - 3

    lines, columns = counts.shape[:2]

    offset = 0

    for i, prt_val in enumerate(prt):
        if prt_val < 50:
            offset = i
            break

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

    tprt = (cal.d[iprt, 0] + prt *
            (cal.d[iprt, 1] + prt *
             (cal.d[iprt, 2] + prt *
              (cal.d[iprt, 3] + prt *
               (cal.d[iprt, 4])))))

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

    # calibrating thermal channel

    tBB = new_tprt
    tsBB = cal.a[chan] + cal.b[chan] * tBB
    nBB_num = (1.1910427 * 0.000010) * cal.c_wn[chan] ** 3
    nBB = nBB_num / (np.exp((1.4387752 * cal.c_wn[chan]) / tsBB) - 1.0)

    Nlin = (cal.n_s[chan] +
            (((nBB - cal.n_s[chan])
              * (new_space - counts.astype(float)))
             / (new_space - new_ict)))
    Ncor = cal.b0[chan] + Nlin * (cal.b1[chan] + cal.b2[chan] * Nlin)
    Ne = Ncor
    tsE = ((1.4387752 * cal.c_wn[chan])
           / np.log(1.0 + nBB_num / Ne))
    bt = (tsE - cal.a[chan]) / cal.b[chan]

    if chan == 0:
        bt = np.where((counts - new_space) >= 0, 0.0, bt)

    # Mask values outside valid range
    bt = np.where(np.logical_or(bt < 170.0, bt > 350.0), np.nan, bt)

    return bt
