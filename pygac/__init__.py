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

from pygac.version import __version__  # noqa
from pygac.reader import ReaderError
from pygac.gac_klm import GACKLMReader
from pygac.gac_pod import GACPODReader
from pygac.lac_klm import LACKLMReader
from pygac.lac_pod import LACPODReader
from pygac.utils import file_opener

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


_Readers = [GACKLMReader, LACKLMReader, GACPODReader, LACPODReader]

def get_reader(filename, fileobj=None):
    """Read the GAC/LAC KLM/POD data.

        Args:
            filename (str): Path to GAC/LAC file
            fileobj: An open file object to read from. (optional)
        """
    found_reader = False
    with file_opener(fileobj or filename) as open_file:
        for i, Reader in enumerate(_Readers):
            try:
                reader = Reader.fromfile(filename, fileobj=open_file)
                found_reader = True
                index = i
            except ReaderError:
                LOG.debug("%s failed to read the file!" % Reader.__name__)
            finally:
                open_file.seek(0)
            if found_reader:
                break
    if not found_reader:
        raise ValueError('Unable to read the file "%s"' % filename)
    # Move the Reader in front of _Readers. Chance is high that the 
    # next file is of the same kind.
    _Readers.insert(0, _Readers.pop(index))
    return reader


def l1c_processor(filename, start_line, end_line, fileobj=None):
    tic = datetime.datetime.now()
    LOG.info("Building file: %s", str(filename))
    reader = get_reader(filename, fileobj=fileobj)
    reader.save(start_line, end_line)
    LOG.info("Processing took: %s", str(datetime.datetime.now() - tic))

