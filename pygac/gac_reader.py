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

from pygac.reader import Reader, ReaderError
import pygac.pygac_geotiepoints as gtp


LOG = logging.getLogger(__name__)


class GACReader(Reader):
    """Reader for GAC data."""

    # Scanning frequency (scanlines per millisecond)
    scan_freq = 2.0 / 1000.0
    # Max scanlines
    max_scanlines = 15000

    def __init__(self, *args, **kwargs):
        """Init the GAC reader."""
        super(GACReader, self).__init__(*args, **kwargs)
        self.scan_width = 409
        self.lonlat_interpolator = gtp.gac_lat_lon_interpolator

    @classmethod
    def _validate_header(cls, header):
        """Check if the header belongs to this reader."""
        # call super to enter the Method Resolution Order (MRO)
        super(GACReader, cls)._validate_header(header)
        LOG.debug("validate header")
        data_set_name = header['data_set_name'].decode()
        # split header into parts
        creation_site, transfer_mode, platform_id = (
            data_set_name.split('.')[:3])
        if transfer_mode != 'GHRR':
            raise ReaderError('Improper transfer mode "%s"!' % transfer_mode)
