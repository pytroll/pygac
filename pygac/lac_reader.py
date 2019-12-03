#!/usr/bin/python
# Copyright (c) 2014-2019.
#

# Author(s):

#   Abhay Devasthale <abhay.devasthale@smhi.se>
#   Adam Dybbroe <adam.dybbroe@smhi.se>
#   Sajid Pareeth <sajid.pareeth@fmach.it>
#   Martin Raspaud <martin.raspaud@smhi.se>

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
"""The LAC reader."""

from pygac.reader import Reader
import pygac.pygac_geotiepoints as gtp


class LACReader(Reader):
    """Reader for LAC data."""

    # Scanning frequency (scanlines per millisecond)
    scan_freq = 6.0 / 1000.0

    def __init__(self, *args, **kwargs):
        """Init the LAC reader."""
        super(LACReader, self).__init__(*args, **kwargs)
        self.scan_width = 2048
        self.lonlat_interpolator = gtp.lac_lat_lon_interpolator
