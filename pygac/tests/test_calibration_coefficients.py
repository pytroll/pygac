#!/usr/bin/env python

# Copyright (c) 2014-2020 Pytroll Developers

# Author(s):

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

"""Test function for the calibration coeffictions handling.
"""

import sys
import unittest
try:
    import mock
except ImportError:
    from unittest import mock
import numpy as np

from pygac.calibration import calibrate_solar, Calibrator
from pygac.configuration import _config

# dummy user json file including only noaa19 data with changed channel 1 coefficients
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
class TestCalibrationCoefficientsHandling(unittest.TestCase):

    @mock.patch('pygac.calibration.open', mock.mock_open(read_data=user_json_file))
    def test_user_coefficients_file(self):
        # reset coefficients
        Calibrator.coeffs = None
        Calibrator._version = None

        if sys.version_info.major < 3:
            cal = Calibrator('noaa19')
        else:
            with self.assertWarnsRegex(RuntimeWarning,
                                       "Unknown calibration coefficients version!"):
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

        # re-reset coefficients
        Calibrator.coeffs = None
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
    mysuite.addTest(loader.loadTestsFromTestCase(TestCalibrationCoefficientsHandling))

    return mysuite


if __name__ == '__main__':
    unittest.main()
