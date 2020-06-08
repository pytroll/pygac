#!/usr/bin/env python

# Copyright (c) 2014-2019 Pytroll Developers

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

"""Test function for the POD calibration.
"""

import sys
import unittest
try:
    import mock
except ImportError:
    from unittest import mock
import numpy as np

from pygac.calibration import calibrate_solar, calibrate_thermal, Calibrator
from pygac.configuration import _config

# user json file including only noaa19 data with changed channel 1 coefficients
user_json_file = """{
    "noaa19": {
        "channel_1": {
            "dark_count": 0,
            "gain_switch": 1000,
            "s0": 2,
            "s1": 0,
            "s2": 0
        },
        "channel_2": {
            "dark_count": 39.0,
            "gain_switch": 500.37,
            "s0": 0.122,
            "s1": 0.95,
            "s2": -0.039
        },
        "channel_3a": {
            "dark_count": 39.4,
            "gain_switch": 496.11,
            "s0": 0.1,
            "s1": 0.0,
            "s2": 0.0
        },
        "channel_3b": {
            "b0": 0.0,
            "b1": 0.0,
            "b2": 0.0,
            "centroid_wavenumber": 2670.2425,
            "space_radiance": 0.0,
            "to_eff_blackbody_intercept": 1.6863857,
            "to_eff_blackbody_slope": 0.9974112191806167
        },
        "channel_4": {
            "b0": 5.7,
            "b1": -0.11187000000000002,
            "b2": 0.00054668,
            "centroid_wavenumber": 927.92374,
            "space_radiance": -5.49,
            "to_eff_blackbody_intercept": 0.39419031,
            "to_eff_blackbody_slope": 0.9986718662850276
        },
        "channel_5": {
            "b0": 3.58,
            "b1": -0.05991000000000002,
            "b2": 0.00024985,
            "centroid_wavenumber": 831.28619,
            "space_radiance": -3.39,
            "to_eff_blackbody_intercept": 0.2636462,
            "to_eff_blackbody_slope": 0.9990463103920997
        },
        "date_of_launch": "2009-02-05T00:57:36.000000Z",
        "thermometer_1": {
            "d0": 276.6067,
            "d1": 0.051111,
            "d2": 1.405783e-06,
            "d3": 0.0,
            "d4": 0.0
        },
        "thermometer_2": {
            "d0": 276.6119,
            "d1": 0.05109,
            "d2": 1.496037e-06,
            "d3": 0.0,
            "d4": 0.0
        },
        "thermometer_3": {
            "d0": 276.6311,
            "d1": 0.051033,
            "d2": 1.49699e-06,
            "d3": 0.0,
            "d4": 0.0
        },
        "thermometer_4": {
            "d0": 276.6268,
            "d1": 0.051058,
            "d2": 1.49311e-06,
            "d3": 0.0,
            "d4": 0.0
        }
    }
}"""


@mock.patch('pygac.configuration.get_config', mock.Mock(return_value=_config))
class TestGenericCalibration(unittest.TestCase):

    def test_calibration_vis(self):
        counts = np.array([[0, 0, 0, 0, 0,
                            512, 512, 512, 512, 512,
                            1023, 1023, 1023, 1023, 1023],
                           [41, 41, 41, 41, 41,
                            150, 150, 150, 150, 150,
                            700, 700, 700, 700, 700]])
        year = 2010
        jday = 1
        spacecraft_id = "noaa19"
        corr = 1
        channel = 0

        ref1 = calibrate_solar(counts[:, channel::5], channel, year, jday,
                               spacecraft_id, corr)

        channel = 1

        ref2 = calibrate_solar(counts[:, channel::5], channel, year, jday,
                               spacecraft_id, corr)

        channel = 2

        data = np.ma.array(counts[:, channel::5], mask=True)

        ref3 = calibrate_solar(data, channel, year, jday,
                               spacecraft_id, corr)

        expected = (np.array([[np.nan, 27.37909518, 110.60103456],
                              [0.11943135, 6.03671211, 57.99695154]]),
                    np.array([[np.nan, 3.05229160e+01, 1.24811455e+02],
                              [1.23011792e-01, 6.82715447e+00, 6.52122414e+01]]),
                    np.array([[0., 523.41775, 1034.41775],
                              [41., 150., 711.41775]]))
        np.testing.assert_allclose(ref1, expected[0])
        np.testing.assert_allclose(ref2, expected[1])
        np.testing.assert_allclose(ref3, expected[2])

    def test_calibration_ir(self):
        counts = np.array([[0, 0, 612, 0, 0,
                            512, 512, 487, 512, 512,
                            923, 923, 687, 923, 923],
                           [41, 41, 634, 41, 41,
                            150, 150, 461, 150, 150,
                            700, 700, 670, 700, 700],
                           [241, 241, 656, 241, 241,
                            350, 350, 490, 350, 350,
                            600, 600, 475, 600, 600]])
        prt_counts = np.array([0, 230, 230])
        ict_counts = np.array([[745.3, 397.9, 377.8],
                               [744.8, 398.1, 378.4],
                               [745.7, 398., 378.3]])
        space_counts = np.array([[987.3,  992.5,  989.4],
                                 [986.9,  992.8,  989.6],
                                 [986.3,  992.3,  988.9]])

        spacecraft_id = "noaa19"
        ch3 = calibrate_thermal(counts[:, 2::5],
                                prt_counts,
                                ict_counts[:, 0],
                                space_counts[:, 0],
                                line_numbers=np.array([1, 2, 3]),
                                channel=3,
                                spacecraft=spacecraft_id)

        expected_ch3 = np.array([[298.36772477, 305.24899954, 293.23847375],
                                 [296.96053595, 306.49432811, 294.48914038],
                                 [295.47715016, 305.10182601, 305.83036782]])

        np.testing.assert_allclose(expected_ch3, ch3)

        ch4 = calibrate_thermal(counts[:, 3::5],
                                prt_counts,
                                ict_counts[:, 1],
                                space_counts[:, 1],
                                line_numbers=np.array([1, 2, 3]),
                                channel=4,
                                spacecraft=spacecraft_id)

        expected_ch4 = np.array([[326.57669548, 275.34893211, 197.68844955],
                                 [323.01324859, 313.20717645, 249.3633716],
                                 [304.58097221, 293.57932356, 264.0630027]])

        np.testing.assert_allclose(expected_ch4, ch4)

        ch5 = calibrate_thermal(counts[:, 4::5],
                                prt_counts,
                                ict_counts[:, 2],
                                space_counts[:, 2],
                                line_numbers=np.array([1, 2, 3]),
                                channel=5,
                                spacecraft=spacecraft_id)

        expected_ch5 = np.array([[326.96168274, 272.09013413, 188.26784127],
                                 [323.15638147, 312.67331324, 244.18437795],
                                 [303.43940924, 291.64944851, 259.97304154]])

        np.testing.assert_allclose(expected_ch5, ch5)

    @mock.patch('pygac.calibration.open', mock.mock_open(read_data=user_json_file))
    def test_user_coefficients_file(self):
        # reset coefficients
        Calibrator.default_coeffs = None
        Calibrator._version = None

        if sys.version_info.major < 3:
            cal = Calibrator('noaa19')
        else:
            with self.assertWarnsRegex(RuntimeWarning,
                                       "Unknown default calibration coefficients version!"):
                cal = Calibrator('noaa19')

        self.assertEqual(cal.dark_count[0], 0)
        self.assertEqual(cal.gain_switch[0], 1000)
        self.assertEqual(cal.s0[0], 2)
        self.assertEqual(cal.s1[0], 0)
        self.assertEqual(cal.s2[0], 0)
        # The coefficients are choosen to preserve the counts of channel 1 if counts are less than 1000
        counts = np.arange(10)
        year = 2010
        jday = 1
        spacecraft_id = "noaa19"
        channel = 0
        scaled_radiance = calibrate_solar(counts, channel, year, jday, spacecraft_id)
        np.testing.assert_allclose(scaled_radiance, counts)

        # reset coefficients
        Calibrator.default_coeffs = None
        Calibrator._version = None

    def test_custom_coefficients(self):
        custom_coeffs = {
            "channel_1": {
                "dark_count": 0,
                "gain_switch": 1000,
                "s0": 2,
                "s1": 0,
                "s2": 0
            }
        }
        # The coefficients are choosen to preserve the counts of channel 1 if counts are less than 1000
        counts = np.arange(10)
        year = 2010
        jday = 1
        spacecraft_id = "noaa19"
        channel = 0
        scaled_radiance = calibrate_solar(counts, channel, year, jday, spacecraft_id,
                                          custom_coeffs=custom_coeffs)
        np.testing.assert_allclose(scaled_radiance, counts)

    @unittest.skipIf(sys.version_info.major < 3, "Skipped in python2!")
    def test_vis_deprecation_warning(self):
        counts = np.arange(10)
        year = 2010
        jday = 1
        spacecraft_id = "noaa19"
        channel = 0
        corr = 2
        message = (
            "Using the 'corr' argument is depricated in favor of making the units"
            " of the function result clear. Please make any unit conversion outside this function."
        )
        with self.assertWarnsRegex(DeprecationWarning, message):
            calibrate_solar(counts, channel, year, jday, spacecraft_id, corr=corr)

def suite():
    """The suite for test_slerp
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestGenericCalibration))

    return mysuite


if __name__ == '__main__':
    unittest.main()
