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

from pygac.version import __version__  # noqa

LOG = logging.getLogger(__name__)
try:
    CONFIG_FILE = os.environ['PYGAC_CONFIG_FILE']
except KeyError:
    LOG.error('Environment variable PYGAC_CONFIG_FILE not set!')
    CONFIG_FILE = ''

if not os.path.exists(CONFIG_FILE) or not os.path.isfile(CONFIG_FILE):
    LOG.warning(
        str(CONFIG_FILE) + " pointed to by the environment "
        + "variable PYGAC_CONFIG_FILE is not a file or does not exist!")


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
