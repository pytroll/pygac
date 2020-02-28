#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c14526.ad.smhi.se>

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

import logging
import os
import numpy as np
try:
    import ConfigParser
except ImportError:
    import configparser as ConfigParser

from pygac.version import __version__  # noqa

LOG = logging.getLogger(__name__)

_CONFIG_FILE = ''


def get_config_file():
    """Return the config file path."""
    global _CONFIG_FILE
    if not _CONFIG_FILE:
        try:
            LOG.info('Config file was not explicitly set. Use'
                     ' environment variable "PYGAC_CONFIG_FILE"')
            _CONFIG_FILE = os.environ["PYGAC_CONFIG_FILE"]
        except KeyError:
            LOG.error('Environment variable PYGAC_CONFIG_FILE not set!')
    if not os.path.isfile(_CONFIG_FILE):
        raise FileNotFoundError('Given config path "%s" is not a file!'
                                % _CONFIG_FILE)
    return _CONFIG_FILE


def set_config_file(path):
    """Set the module config file."""
    global _CONFIG_FILE
    _CONFIG_FILE = str(path)


def get_config():
    """Retrun the module configuration."""
    config_file = get_config_file()
    config = ConfigParser.ConfigParser()
    try:
        config.read(config_file)
    except ConfigParser.NoSectionError as exception:
        LOG.error('Failed reading configuration file: "%s"' % config_file)
        raise exception
    return config


def get_absolute_azimuth_angle_diff(sat_azi, sun_azi):
    """Calculates absolute azimuth difference angle. """
    rel_azi = abs(sat_azi - sun_azi)
    rel_azi = rel_azi % 360
    # Not using np.where to avoid copying array
    rel_azi[rel_azi > 180] = 360.0 - rel_azi[rel_azi > 180]
    return rel_azi


def centered_modulus(array, divisor):
    """Transform array to half open range ]-divisor/2, divisor/2]."""
    arr = array % divisor
    arr[arr > divisor / 2] -= divisor
    return arr


def calculate_sun_earth_distance_correction(jday):
    """Calculate the sun earth distance correction.

    In 2008 3-4 different equations of ESD were considered.
    This one was chosen as it at the time gave reflectances most closely
    matching the PATMOS-x data provided then by Andy Heidinger.

    Formula might need to be reconsidered if jday is updated to a float.

    """
    # Earth-Sun distance correction factor
    corr = 1.0 - 0.0334 * np.cos(2.0 * np.pi * (jday - 2) / 365.25)
    return corr
