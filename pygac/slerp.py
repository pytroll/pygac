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

"""slerp implementation in numpy
"""


import unittest
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
    omega = np.arccos(dot(cp0, cp1))[:, :, np.newaxis]
    print omega.shape, t.shape
    return toll(np.sin((1 - t) * omega) / np.sin(omega) * cp0 +
                np.sin(t * omega) / np.sin(omega) * cp1)


class TestSlerp(unittest.TestCase):

    def test_slerp(self):
        lon0, lat0 = (0, 0)
        lon1, lat1 = (0, 1)
        self.assertTrue(
            np.allclose(slerp(lon0, lat0, lon1, lat1, 0.5), (0, 0.5)))

    def test_slerp_datum(self):
        lon0, lat0 = (183, 0)
        lon1, lat1 = (179, 0)
        res = slerp(lon0, lat0, lon1, lat1, 0.5)
        res %= 360
        self.assertTrue(
            np.allclose(res, (181, 0)))

    def test_slerp_pole(self):
        lon0, lat0 = (0, 89)
        lon1, lat1 = (180, 89)
        res = slerp(lon0, lat0, lon1, lat1, 0.5)
        self.assertTrue(
            np.allclose(res[:, :, 1], 90))

        lon0, lat0 = (-90, 89)
        lon1, lat1 = (90, 89)
        res = slerp(lon0, lat0, lon1, lat1, 0.5)
        self.assertTrue(
            np.allclose(res[:, :, 1], 90))

        lon0, lat0 = (0, 89)
        lon1, lat1 = (180, 87)
        res = slerp(lon0, lat0, lon1, lat1, 0.5)
        self.assertTrue(
            np.allclose(res, (180, 89)))

    def test_slerp_vec(self):
        lon0 = np.array([[0, 0],
                         [0, 0]])
        lat0 = np.array([[0, 0],
                         [0, 0]])
        lon1 = np.array([[0, 0],
                         [0, 0]])
        lat1 = np.array([[1, 1],
                         [1, 1]])

        res = slerp(lon0, lat0, lon1, lat1, 0.5)
        self.assertTrue(np.allclose(res[:, :, 0], 0))
        self.assertTrue(np.allclose(res[:, :, 1], 0.5))

        lon0 = np.array([[183, 0],
                         [-90, 0]])
        lat0 = np.array([[0, 89],
                         [89, 89]])
        lon1 = np.array([[179, 180],
                         [90, 180]])
        lat1 = np.array([[0, 89],
                         [89, 87]])

        res = slerp(lon0, lat0, lon1, lat1, 0.5)

        self.assertTrue(np.allclose(res[0, 0, :] % 360, (181, 0)))
        self.assertTrue(np.allclose(res[0, 1, 1], 90))
        self.assertTrue(np.allclose(res[1, 0, 1], 90))
        self.assertTrue(np.allclose(res[1, 1, :], (180, 89)))

    def test_slerp_tvec(self):
        lon0 = np.array([[0, 0],
                         [0, 0],
                         [0, 0],
                         [0, 0],
                         [0, 0],
                         [0, 0],
                         [0, 0]])
        lat0 = np.array([[0, 0],
                         [5, 0],
                         [10, 0],
                         [15, 0],
                         [20, 0],
                         [25, 0],
                         [30, 0]])
        lon1 = np.array([[0, 0],
                         [0, 0],
                         [0, 0],
                         [0, 0],
                         [0, 0],
                         [0, 0],
                         [0, 0]])
        lat1 = np.array([[45, 45],
                         [45, 45],
                         [45, 45],
                         [45, 45],
                         [45, 45],
                         [45, 45],
                         [45, 45]])

        t = np.array([[0.5, 0, 0.2, 0.4, 0.6, 0.8, 1]]).T
        t = t[:, :, np.newaxis]
        res = slerp(lon0, lat0, lon1, lat1, t)
        expected = np.array([[22.5, 22.5],
                             [5., 0.],
                             [17., 9.],
                             [27., 18.],
                             [35., 27.],
                             [41., 36.],
                             [45., 45.]])

        self.assertTrue(np.allclose(res[:, :, 0], 0))
        self.assertTrue(np.allclose(res[:, :, 1], expected))

    def test_slerp_interp_cols(self):
        lons = np.array([[2, 4, 6],
                         [1, 3, 5]])
        lats = np.array([[0, 0, 0],
                         [0, 0, 0]])
        cols_subset = np.array([2, 4, 6])
        cols_full = np.arange(8)

        res = slerp_interp_cols(lons, lats, cols_subset, cols_full)
        # print res

if __name__ == '__main__':
    unittest.main()
