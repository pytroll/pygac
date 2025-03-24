#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014, 2015 Martin Raspaud

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
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

"""Generic reader for GAC data.

Can't be used as is, has to be subclassed to add specific read functions.
"""

import logging
import warnings

try:
    from pyorbital.geoloc_instrument_definitions import avhrr_gac_from_times
except ImportError:
    # In pyorbital 1.9.2 and earlier avhrr_gac returned LAC geometry
    from pyorbital.geoloc_instrument_definitions import avhrr_gac
    def avhrr_gac_from_times(times, points):
        return avhrr_gac(times, points*5+3.5)
    warnings.warn('pyorbital version does not support avhrr_gac_from_times. ' +
                  'Computation of missing longitude/latitudes may be incorrect.')


import pygac.pygac_geotiepoints as gtp
from pygac.pygac_geotiepoints import GAC_LONLAT_SAMPLE_POINTS
from pygac.reader import Reader, ReaderError

LOG = logging.getLogger(__name__)


class GACReader(Reader):
    """Reader for GAC data."""

    # Scanning frequency (scanlines per millisecond)
    scan_freq = 2.0 / 1000.0
    # Max scanlines
    max_scanlines = 15000
    lonlat_sample_points = GAC_LONLAT_SAMPLE_POINTS

    def __init__(self, *args, **kwargs):
        """Init the GAC reader."""
        super().__init__(*args, **kwargs)
        self.scan_width = 409
        self.lonlat_interpolator = gtp.gac_lat_lon_interpolator
        self.geoloc_definition = avhrr_gac_from_times

    @classmethod
    def _validate_header(cls, header):
        """Check if the header belongs to this reader."""
        # call super to enter the Method Resolution Order (MRO)
        super(GACReader, cls)._validate_header(header)
        LOG.debug("validate header")
        data_set_name = header["data_set_name"].decode()
        # split header into parts
        creation_site, transfer_mode, platform_id = (
            data_set_name.split(".")[:3])
        if transfer_mode != "GHRR":
            raise ReaderError('Improper transfer mode "%s"!' % transfer_mode)
