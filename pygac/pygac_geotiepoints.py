#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2010-2012.

# Author(s):

#   Adam Dybbroe <adam.dybbroe@smhise>
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

"""Interpolation of geographical tiepoints."""


import geotiepoints as gtp
import numpy as np

GAC_LONLAT_SAMPLE_POINTS = np.arange(4, 405, 8)
LAC_LONLAT_SAMPLE_POINTS = np.arange(24, 2048, 40)


def gac_lat_lon_interpolator(lons_subset, lats_subset):
    """Interpolate lat-lon values in the AVHRR GAC data.

    Each GAC row has 409 pixels.
    But lat-lon values are provided for every eigth pixel,
    ranging from 5 to 405. Interpolate to full resolution.

    """
    cols_full = np.arange(409)
    return lat_lon_interpolator(lons_subset, lats_subset, GAC_LONLAT_SAMPLE_POINTS, cols_full)



def lac_lat_lon_interpolator(lons_subset, lats_subset):
    """Interpolate lat-lon values in the AVHRR LAC data."""
    cols_full = np.arange(2048)
    return lat_lon_interpolator(lons_subset, lats_subset, LAC_LONLAT_SAMPLE_POINTS, cols_full)


def lat_lon_interpolator(lons_subset, lats_subset, cols_subset, cols_full):
    """Interpolate lat-lon values in the AVHRR data."""
    lines = lats_subset.shape[0]
    rows_subset = np.arange(lines)
    rows_full = np.arange(lines)

    along_track_order = 1
    cross_track_order = 3

    satint = gtp.SatelliteInterpolator((lons_subset, lats_subset),
                                       (rows_subset, cols_subset),
                                       (rows_full, cols_full),
                                       along_track_order,
                                       cross_track_order)

    return satint.interpolate()
