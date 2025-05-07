#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 Martin Raspaud

# Author(s):

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

"""slerp implementation in numpy"""

import numpy as np


def tocart(lon, lat):
    rlat = np.deg2rad(lat)
    rlon = np.deg2rad(lon)
    x__ = np.cos(rlat) * np.cos(rlon)
    y__ = np.cos(rlat) * np.sin(rlon)
    z__ = np.sin(rlat)
    return np.dstack((x__, y__, z__))


def toll(arr):
    lat = np.arcsin(arr.take(2, -1))
    lon = np.arctan2(arr.take(1, -1), arr.take(0, -1))
    return np.rad2deg(np.dstack((lon.squeeze(), lat.squeeze())))


def dot(a, b):
    return np.sum(a * b, -1)


def slerp(lon0, lat0, lon1, lat1, t):
    cp0 = tocart(lon0, lat0)
    cp1 = tocart(lon1, lat1)

    dot_product = np.clip(dot(cp0, cp1), -1.0, 1.0)
    identical_mask = np.isclose(dot_product, 1.0)

    omega = np.arccos(dot_product)[:, :, np.newaxis]
    sin_omega = np.sin(omega)

    interp = (
        np.sin((1 - t) * omega) / np.where(sin_omega == 0, 1, sin_omega) * cp0
        + np.sin(t * omega) / np.where(sin_omega == 0, 1, sin_omega) * cp1
    )

    interp[identical_mask] = cp0[identical_mask]
    return toll(interp)
