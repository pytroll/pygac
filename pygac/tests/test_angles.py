"""Test function for the angle calculation."""

import unittest

import numpy as np

from pygac.utils import centered_modulus, get_absolute_azimuth_angle_diff


class TestAngles(unittest.TestCase):
    """Test function for the angle calculation."""

    def test_azidiff_angle(self):
        """Test function for the azidiff angle."""
        sat_az = np.ma.array([[48.0,  56.0, 64.0,  72.0],
                              [80.0,  88.0, 96.0, 104.0],
                              [-80.0,  -88.0, -96.0, -104.0],
                              [-180.0,  -188.0, -196.0, -204.0]], mask=False)
        sun_az = np.ma.array([[148.0,  156.0, 164.0,  172.0],
                              [180.0,  188.0, 196.0, 204.0],
                              [180.0,  188.0, 196.0, 204.0],
                              [185.0,  193.0, 201.0, 209.0]], mask=False)

        res = np.ma.array([[100., 100., 100., 100.],
                           [100., 100., 100., 100.],
                           [100.,  84.,  68.,  52.],
                           [5.,  21.,  37.,  53.]],
                          mask=False)
        rel_azi = get_absolute_azimuth_angle_diff(sat_az, sun_az)

        np.testing.assert_allclose(rel_azi, res)

    def test_centered_modulus(self):
        """Test centered_modulus."""
        angles = np.ma.array(
            [[180.0,  -180.0, -179.9,  201.0],
             [80.0,  360.0, -360.0, 604.0],
             [-80.0,  -88.0, -796.0, -104.0],
             [-3.0,  -188.0, -196.0, -204.0]], mask=False)
        expected = np.ma.array(
            [[180.0,  180.0, -179.9,  -159.0],
             [80.0,  0.0, 0.0, -116.0],
             [-80.0,  -88.0, -76.0, -104.0],
             [-3.0,  172.0, 164.0, 156.0]], mask=False)
        transformed = centered_modulus(angles, 360.0)
        np.testing.assert_allclose(transformed, expected)
