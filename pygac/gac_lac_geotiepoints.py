#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2010-2015.

# Author(s):

#   Abhay Devasthale <abhay.devasthale@smhi.se>
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

"""Interpolation of geographical tiepoints for GAC and LAC data.
"""

import numpy as np

from geotiepoints import SatelliteInterpolator


def geo_interpolator(lons_subset, lats_subset, scan_width, scan_points):
    """ interpolates lat-lon values in the AVHRR GAC data from every eigth pixel to all pixels
        Each GAC row has total 409 pixels. But lat-lon values are provided for every eigth pixel starting from pixel 5 and ending at pixel 405 """

    cols_subset = scan_points
    cols_full = np.arange(scan_width)
    lines = lats_subset.shape[0]
    rows_subset = np.arange(lines)
    rows_full = np.arange(lines)

    along_track_order = 1
    cross_track_order = 3

    satint = SatelliteInterpolator((lons_subset, lats_subset),
                                   (rows_subset, cols_subset),
                                   (rows_full, cols_full),
                                   along_track_order,
                                   cross_track_order)

    return satint.interpolate()
