"""Test the slerp module
"""


import unittest

import numpy as np

from pygac.slerp import slerp


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
