#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014-2015, 2019 Pytroll Developers

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
#   Abhay Devasthale <abhay.devasthale@smhi.se>

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

"""Calibration coefficients and generic calibration functions
"""
from __future__ import division
import numpy as np

coeffs = {
    'metopb': {'ah': np.array([0.166, 0.183, 0.201]),
               'al': np.array([0.055, 0.061, 0.029]),
               'bh': np.array([2.019, 1.476, 1.478]),
               'bl': np.array([2.019, 1.476, 1.478]),
               'ch': np.array([-0.201, -0.137, -0.033]),
               'cl': np.array([-0.201, -0.137, -0.033]),
               'c_dark': np.array([39.70, 40.00, 40.30]),
               'c_s': np.array([501.12, 500.82, 501.32]),
               'l_date': 2012.77,
               'd': np.array([[0, 0, 0, 0, 0],  # reset prt
                              [276.6194, 0.050919, 1.471E-06, 0.0, 0.0],
                              [276.6511, 0.050892, 1.489E-06, 0.0, 0.0],
                              [276.6597, 0.050845, 1.521E-06, 0.0, 0.0],
                              [276.3685, 0.050992, 1.482E-06, 0.0, 0.0]]),
               'n_s': np.array([0.0, -4.98, -3.40]),
               'c_wn': np.array([2664.3384, 933.71521, 839.72764]),
               'a': np.array([1.7711318, 0.51860807, 0.40059787]),
               'b': np.array([1.0 / 1.0029931,
                              1.0 / 1.0013778,
                              1.0 / 1.0011702]),
               'b0': np.array([0.0, 5.44, 3.84]),
               'b1': np.array([1 - 0.0, 1 - 0.10152, 1 - 0.06249]),
               'b2': np.array([0.0, 0.00046964, 0.00025239]),
               },
    'metopa': {'ah': np.array([0.169, 0.199, 0.213]),
               'al': np.array([0.056, 0.066, 0.030]),
               'bh': np.array([0.609, 0.980, -0.016]),
               'bl': np.array([0.609, 0.980, -0.016]),
               'ch': np.array([-0.029, -0.016, -0.033]),
               'cl': np.array([-0.029, -0.016, -0.033]),
               'c_dark': np.array([40.43, 39.75, 41.8]),
               'c_s': np.array([501.0, 500.0, 502.0]),
               'l_date': 2006.7995,
               'd': np.array([[0, 0, 0, 0, 0],  # reset prt
                              [276.6194, 0.050919, 1.471E-06, 0.0, 0.0],
                              [276.6511, 0.050892, 1.489E-06, 0.0, 0.0],
                              [276.6597, 0.050845, 1.521E-06, 0.0, 0.0],
                              [276.3685, 0.050992, 1.482E-06, 0.0, 0.0]]),
               'n_s': np.array([0.0, -4.98, -3.40]),
               'c_wn': np.array([2687.0392, 927.27630, 837.80762]),
               'a': np.array([2.0653147, 0.56503332, 0.38472766]),
               'b': np.array([1.0 / 1.0034418,
                              1.0 / 1.0015090,
                              1.0 / 1.0011264]),
               'b0': np.array([0.0, 5.44, 3.84]),
               'b1': np.array([1 - 0.0, 1 - 0.10152, 1 - 0.06249]),
               'b2': np.array([0.0, 0.00046964, 0.00025239]),
               },
    'tirosn': {'ah': np.array([0.115, 0.133, 0.1]),
               'al': np.array([0.115, 0.133, 0.1]),
               'bh': np.array([5.110, 0.717, 0.0]),
               'bl': np.array([5.110, 0.717, 0.0]),
               'ch': np.array([0.0, 0.0, 0.0]),
               'cl': np.array([0.0, 0.0, 0.0]),
               'c_dark': np.array([39.44, 39.40, 37.51]),
               'l_date': 1978.783,
               'd': np.array([[0, 0, 0, 0, 0],  # reset prt
                              [277.659, 0.051275, 1.363e-06, 0, 0],
                              [276.659, 0.051275, 1.363e-06, 0, 0],
                              [276.659, 0.051275, 1.363e-06, 0, 0],
                              [276.659, 0.051275, 1.363e-06, 0, 0]]),
               'n_s': np.array([-0.0039, -8.130, -8.130]),
               'c_wn': np.array([2655.7409, 913.0537, 913.0537]),
               'a': np.array([1.6485446, 0.53135445, 0.53135445]),
               'b': np.array([1.0 / 1.0020894,
                              1.0 / 1.0014343,
                              1.0 / 1.0014343]),
               'b1': np.array([1.0 - 0.015, 1.0 - 0.131942, 1.0 - 0.131942]),
               'b2': np.array([0.011, 0.000673193, 0.000673193]),
               'b0': np.array([0.00195, 6.13, 6.13]),
               },
    'noaa6': {'ah': np.array([0.133, 0.128, 0.10]),
              'al': np.array([0.133, 0.128, 0.10]),
              'bh': np.array([0.900, 0.699, 0.0]),
              'bl': np.array([0.900, 0.699, 0.0]),
              'ch': np.array([0.0, 0.0, 0.0]),
              'cl': np.array([0.0, 0.0, 0.0]),
              'c_dark': np.array([39.44, 39.40, 37.51]),
              'l_date': 1979.490,
              'd': np.array([[0, 0, 0, 0, 0],  # reset prt
                             [277.659, 0.051275, 1.363e-06, 0, 0],
                             [276.659, 0.051275, 1.363e-06, 0, 0],
                             [276.659, 0.051275, 1.363e-06, 0, 0],
                             [276.659, 0.051275, 1.363e-06, 0, 0]]),
              'n_s': np.array([0.0, -3.26, -3.26]),
              'c_wn': np.array([2671.5433, 913.46088, 913.46088]),
              'a': np.array([1.76671100, 0.50395970, 0.50395970]),
              'b': np.array([1.0 / 1.0024428,
                             1.0 / 1.0013592,
                             1.0 / 1.0013592]),
              'b1': np.array([1.0, 1.0 - 0.03964, 1.0 - 0.03964]),
              'b2': np.array([0.0, 0.00016925, 0.00016925]),
              'b0': np.array([0.0, 2.24, 2.24]),
              },
    'noaa7': {'ah': np.array([0.115, 0.127, 0.10]),
              'al': np.array([0.115, 0.127, 0.10]),
              'bh': np.array([3.792, 2.685, 0.0]),
              'bl': np.array([3.972, 2.685, 0.0]),
              'ch': np.array([-0.269, -0.101, 0.0]),
              'cl': np.array([-0.269, -0.101, 0.0]),
              'c_dark': np.array([36.0, 37.0, 39.0]),
              'l_date': 1981.4764,
              'd': np.array([[0, 0, 0, 0, 0],  # reset prt
                             [277.099, 5.048E-2, 2.823E-6, 0, 0],
                             [276.734, 5.069E-2, 2.493E-6, 0, 0],
                             [276.876, 5.148E-2, 1.040E-6, 0, 0],
                             [276.160, 5.128E-2, 1.414E-6, 0, 0]]),
              'n_s': np.array([0.0, -5.16, -4.28]),
              'c_wn': np.array([2684.5233, 928.23757, 841.52137]),
              'a': np.array([1.94882690, 0.52807997, 0.40557027]),
              'b': np.array([1.0 / 1.0029260,
                             1.0 / 1.0014039,
                             1.0 / 1.0011789]),
              'b1': np.array([1.0, 0.89783, 0.93683]),
              'b2': np.array([0.0, 0.0004819, 0.0002425]),
              'b0': np.array([0.0, 5.25, 3.93]),
              },
    'noaa8': {'ah': np.array([0.119, 0.136, 0.10]),
              'al': np.array([0.119, 0.136, 0.10]),
              'bh': np.array([6.065, 7.248, 0.0]),
              'bl': np.array([6.065, 7.248, 0.0]),
              'ch': np.array([0.0, 0.0, 0.0]),
              'cl': np.array([0.0, 0.0, 0.0]),
              'c_dark': np.array([39.44, 39.40, 37.51]),
              'l_date': 1983.241,
              'd': np.array([[0, 0, 0, 0, 0],  # reset prt
                             [276.659, 0.051275, 1.363e-06, 0, 0],
                             [276.659, 0.051275, 1.363e-06, 0, 0],
                             [276.659, 0.051275, 1.363e-06, 0, 0],
                             [276.659, 0.051275, 1.363e-06, 0, 0]]),
              'n_s': np.array([0.0, -3.26, -3.26]),
              'c_wn': np.array([2651.3776, 915.30330, 915.30330]),
              'a': np.array([1.76641050, 0.50017997, 0.50017997]),
              'b': np.array([1.0 / 1.0024260,
                             1.0 / 1.0013460,
                             1.0 / 1.0013460]),
              'b1': np.array([1.0, 1.0 - 0.03964, 1.0 - 0.03964]),
              'b2': np.array([0.0, 0.00016925, 0.00016925]),
              'b0': np.array([0.0, 2.24, 2.24]),
              },
    'noaa9': {'ah': np.array([0.108, 0.122, 0.10]),
              'al': np.array([0.108, 0.122, 0.10]),
              'bh': np.array([4.255, 0.310, 0.0]),
              'bl': np.array([4.255, 0.310, 0.0]),
              'ch': np.array([0.640, 0.642, 0.0]),
              'cl': np.array([0.640, 0.642, 0.0]),
              'c_dark': np.array([38.0, 40.0, 38.0]),
              'l_date': 1984.9480,
              'd': np.array([[0, 0, 0, 0, 0],  # reset prt
                             [277.018000, 0.051280, 0.0, 0, 0],
                             [276.750000, 0.051280, 0.0, 0, 0],
                             [276.862000, 0.051280, 0.0, 0, 0],
                             [276.546000, 0.051280, 0.0, 0, 0]]),

              'n_s': np.array([0.0, -5.530, -3.06]),
              'c_wn': np.array([2690.0451, 930.50230, 845.75000]),
              'a': np.array([1.8832662, 0.5115335, 0.3882150]),
              'b': np.array([1.0 / 1.0028978,
                             1.0 / 1.0013570,
                             1.0 / 1.0011210]),
              'b1': np.array([1.0, 0.88643, 0.95311]),
              'b2': np.array([0.0, 0.0006033, 0.0002198]),
              'b0': np.array([0.0, 5.24, 2.42]),
              },
    'noaa10': {'ah': np.array([0.111, 0.137, 0.10]),
               'al': np.array([0.111, 0.137, 0.10]),
               'bh': np.array([6.087, 0.119, 0.0]),
               'bl': np.array([6.087, 0.119, 0.0]),
               'ch': np.array([-1.039, 0.123, 0.0]),
               'cl': np.array([-1.039, 0.123, 0.0]),
               'c_dark': np.array([39.44, 39.40, 37.51]),
               'l_date': 1986.712,
               'd': np.array([[0, 0, 0, 0, 0],  # reset prt
                              [276.659, 0.051275, 1.363e-06, 0, 0],
                              [276.659, 0.051275, 1.363e-06, 0, 0],
                              [276.659, 0.051275, 1.363e-06, 0, 0],
                              [276.659, 0.051275, 1.363e-06, 0, 0]]),
               'n_s': np.array([0.0, -7.27, -7.29]),
               'c_wn': np.array([2672.6164, 910.49626, 910.49626]),
               'a': np.array([1.7986926, 0.45707063, 0.45707063]),
               'b': np.array([1.0 / 1.0026326,
                              1.0 / 1.0012272,
                              1.0 / 1.0012272]),
               'b1': np.array([1.0, 1.0 - 0.1157, 1.0 - 0.1157]),
               'b2': np.array([0.0, 0.0005885, 0.0005882]),
               'b0': np.array([0.0, 5.76, 5.76]),
               },
    'noaa11': {'ah': np.array([0.110, 0.118, 0.0]),
               'al': np.array([0.110, 0.118, 0.0]),
               'bh': np.array([0.632, -0.037, 0.0]),
               'bl': np.array([0.632, -0.037, 0.0]),
               'ch': np.array([-0.044, 0.072, 0.0]),
               'cl': np.array([-0.044, 0.072, 0.0]),
               'c_dark': np.array([40.0, 40.0, 40.0]),
               'l_date': 1988.7310,
               'd': np.array([[0, 0, 0, 0, 0],  # reset prt
                              [276.597, 0.051275, 1.363e-06, 0, 0],
                              [276.597, 0.051275, 1.363e-06, 0, 0],
                              [276.597, 0.051275, 1.363e-06, 0, 0],
                              [276.597, 0.051275, 1.363e-06, 0, 0]]),
               'n_s': np.array([0.0, -8.055, -3.51]),
               'c_wn': np.array([2680.05, 927.462, 840.746]),
               'a': np.array([1.738973, 0.321199, 0.048652]),
               'b': np.array([1.0 / 1.003354,
                              1.0 / 1.001213,
                              1.0 / 1.000664]),
               'b1': np.array([1.0, 0.84120, 0.94598]),
               'b2': np.array([0.0, 0.0008739, 0.0002504]),
               'b0': np.array([0.0, 7.21, 2.92]),
               },
    'noaa12': {'ah': np.array([0.121, 0.148, 0.10]),
               'al': np.array([0.121, 0.148, 0.10]),
               'bh': np.array([2.032, 1.323, 0.0]),
               'bl': np.array([2.032, 1.323, 0.0]),
               'ch': np.array([-0.032, -0.008, 0.0]),
               'cl': np.array([-0.032, -0.008, 0.0]),
               'c_dark': np.array([41.0, 40.0, 40.0]),
               'l_date': 1991.3669,
               'd': np.array([[0, 0, 0, 0, 0],  # reset prt
                              [276.597, 0.051275, 1.363e-06, 0, 0],
                              [276.597, 0.051275, 1.363e-06, 0, 0],
                              [276.597, 0.051275, 1.363e-06, 0, 0],
                              [276.597, 0.051275, 1.363e-06, 0, 0]]),
               'n_s': np.array([0.0, -5.510, -2.51]),
               'c_wn': np.array([2651.7708, 922.36261, 838.02678]),
               'a': np.array([1.90527390, 0.63404209, 0.41086587]),
               'b': np.array([1.0 / 1.0030100,
                              1.0 / 1.0017076,
                              1.0 / 1.0012010]),
               'b1': np.array([1.0, 0.88929, 0.96299]),
               'b2': np.array([0.0, 0.0005968, 0.0001775]),
               'b0': np.array([0.0, 5.11, 1.91]),
               },
    'noaa14': {'ah': np.array([0.121, 0.152, 0.10]),
               'al': np.array([0.121, 0.152, 0.10]),
               'bh': np.array([3.555, 0.254, 0.0]),
               'bl': np.array([3.555, 0.254, 0.0]),
               'ch': np.array([-0.339, 0.201, 0.0]),
               'cl': np.array([-0.339, 0.201, 0.0]),
               'c_dark': np.array([41.0, 41.0, 39.0]),
               'l_date': 1994.9966,
               'd': np.array([[0, 0, 0, 0, 0],  # reset prt
                              [276.597, 0.051275, 1.363e-06, 0, 0],
                              [276.597, 0.051275, 1.363e-06, 0, 0],
                              [276.597, 0.051275, 1.363e-06, 0, 0],
                              [276.597, 0.051275, 1.363e-06, 0, 0]]),
               'n_s': np.array([0.0069, -4.05, -2.29]),
               'c_wn': np.array([2654.25, 928.349, 833.040]),
               'a': np.array([1.885330, 0.308384, 0.022171]),
               'b': np.array([1.0 / 1.003839, 1.0 / 1.001443, 1.0 / 1.000538]),
               'b1': np.array([1.00359, 0.92378, 0.96194]),
               'b2': np.array([0.0, 0.0003822, 0.0001742]),
               'b0': np.array([-0.0031, 3.72, 2.00]),
               },
    'noaa15': {'ah': np.array([0.179, 0.206, 0.175]),
               'al': np.array([0.060, 0.069, 0.025]),
               'bh': np.array([-0.069, 0.339, 0.0]),
               'bl': np.array([-0.069, 0.339, 0.0]),
               'ch': np.array([0.002, -0.010, 0.0]),
               'cl': np.array([0.002, -0.010, 0.0]),
               'c_dark': np.array([39.0, 40.0, 39.0]),
               'c_s': np.array([500.0, 500.0, 500.0]),
               'l_date': 1998.3641,
               'd': np.array([[0, 0, 0, 0, 0],  # reset prt
                              [276.60157, 0.051045, 1.36328E-06, 0.0, 0.0],
                              [276.62531, 0.050909, 1.47266E-06, 0.0, 0.0],
                              [276.67413, 0.050907, 1.47656E-06, 0.0, 0.0],
                              [276.59258, 0.050966, 1.47656E-06, 0.0, 0.0]]),
               'n_s': np.array([0.0, -4.50, -3.61]),
               'c_wn': np.array([2695.9743, 925.4075, 839.8979]),
               'a': np.array([1.624481, 0.338243, 0.304856]),
               'b': np.array([1.0 / 1.001989,
                              1.0 / 1.001283,
                              1.0 / 1.000977]),
               'b0': np.array([0.0, 4.76, 3.83]),
               'b1': np.array([1 - 0.0, 1 - 0.0932, 1 - 0.0659]),
               'b2': np.array([0.0, 0.0004524, 0.0002811]),
               },
    'noaa16': {'ah': np.array([0.165, 0.179, 0.187]),
               'al': np.array([0.055, 0.060, 0.027]),
               'bh': np.array([0.839, 0.786, 0.290]),
               'bl': np.array([0.839, 0.786, 0.290]),
               'ch': np.array([-0.051, -0.031, -0.294]),
               'cl': np.array([-0.051, -0.031, -0.294]),
               'c_dark': np.array([39.3, 38.9, 38.4]),
               'c_s': np.array([498.96, 500.17, 499.43]),
               'l_date': 2000.7228,
               'd': np.array([[0, 0, 0, 0, 0],  # reset prt
                              [276.355, 5.562E-02, -1.590E-05,
                                  2.486E-08, -1.199E-11],
                              [276.142, 5.605E-02, -1.707E-05,
                                  2.595E-08, -1.224E-11],
                              [275.996, 5.486E-02, -1.223E-05,
                                  1.862E-08, -0.853E-11],
                              [276.132, 5.494E-02, -1.344E-05,
                                  2.112E-08, -1.001E-11]]),
               'n_s': np.array([0.0, -2.467, -2.009]),
               'c_wn': np.array([2681.2540, 922.34790, 834.61814]),
               'a': np.array([1.6774586, 0.55636216, 0.41430789]),
               'b': np.array([1.0 / 1.0017316,
                              1.0 / 1.0014921,
                              1.0 / 1.0012166]),
               'b0': np.array([0.0, 2.96, 2.25]),
               'b1': np.array([1 - 0.0, 1 - 0.05411, 1 - 0.03665]),
               'b2': np.array([0.0, 0.00024532, 0.00014854]),
               },
    'noaa17': {'ah': np.array([0.172, 0.210, 0.209]),
               'al': np.array([0.057, 0.070, 0.030]),
               'bh': np.array([1.007, 1.474, 2.787]),
               'bl': np.array([1.007, 1.474, 2.787]),
               'ch': np.array([-0.044, -0.118, -0.292]),
               'cl': np.array([-0.044, -0.118, -0.292]),
               'c_dark': np.array([39.99, 39.09, 42.09]),
               'c_s': np.array([501.12, 500.73, 501.37]),
               'l_date': 2002.47912,
               'd': np.array([[0, 0, 0, 0, 0],  # reset prt
                              [276.628, 0.05098, 1.371E-06, 0.0, 0.0],
                              [276.538, 0.05098, 1.371E-06, 0.0, 0.0],
                              [276.761, 0.05097, 1.369E-06, 0.0, 0.0],
                              [276.660, 0.05100, 1.348E-06, 0.0, 0.0]]),
               'n_s': np.array([0.0, -8.55, -3.97]),
               'c_wn': np.array([2669.1414, 928.29959, 840.20289]),
               'a': np.array([1.70002941, 0.56634758, 0.37264803]),
               'b': np.array([1.0 / 1.0026724,
                              1.0 / 1.0015205,
                              1.0 / 1.0010841]),
               'b0': np.array([0.0, 8.22, 4.31]),
               'b1': np.array([1 - 0.0, 1 - 0.15795, 1 - 0.07318]),
               'b2': np.array([0.0, 0.00075579, 0.00030976]),
               },
    'noaa18': {'ah': np.array([0.171, 0.192, 0.175]),
               'al': np.array([0.057, 0.064, 0.025]),
               'bh': np.array([0.603, 0.632, 0.0]),
               'bl': np.array([0.603, 0.632, 0.0]),
               'ch': np.array([0.0, 0.045, 0.0]),
               'cl': np.array([0.0, 0.045, 0.0]),
               'c_dark': np.array([39.44, 39.40, 37.51]),
               'c_s': np.array([500.54, 500.40, 500.56]),
               'l_date': 2005.3833,
               'd': np.array([[0, 0, 0, 0, 0],  # reset prt
                              [276.601, 0.05090, 1.657E-06, 0.0, 0.0],
                              [276.683, 0.05101, 1.482E-06, 0.0, 0.0],
                              [276.565, 0.05117, 1.313E-06, 0.0, 0.0],
                              [276.615, 0.05103, 1.484E-06, 0.0, 0.0]]),
               'n_s': np.array([0.0, -5.53, -2.22]),
               'c_wn': np.array([2660.6468, 928.73452, 834.08306]),
               'a': np.array([1.7222650, 0.54696239, 0.39938376]),
               'b': np.array([1.0 / 1.0028633,
                              1.0 / 1.0014581,
                              1.0 / 1.0011724]),
               'b0': np.array([0.0, 5.82, 2.67]),
               'b1': np.array([1 - 0.0, 1 - 0.11069, 1 - 0.04360]),
               'b2': np.array([0.0, 0.00052337, 0.00017715]),
               },
    'noaa19': {'ah': np.array([0.162, 0.183, 0.175]),
               'al': np.array([0.054, 0.061, 0.025]),
               'bh': np.array([0.626, 0.950, 0.0]),
               'bl': np.array([0.626, 0.950, 0.0]),
               'ch': np.array([-0.044, -0.039, 0.0]),
               'cl': np.array([-0.044, -0.039, 0.0]),
               'c_dark': np.array([38.8, 39.00, 39.4]),
               'c_s': np.array([496.43, 500.37, 496.11]),
               'l_date': 2009.096,
               'd': np.array([[0, 0, 0, 0, 0],  # reset prt
                              [276.6067, 0.051111, 1.405783e-06, 0, 0],
                              [276.6119, 0.051090, 1.496037e-06, 0, 0],
                              [276.6311, 0.051033, 1.496990e-06, 0, 0],
                              [276.6268, 0.051058, 1.493110e-06, 0, 0]]),
               'n_s': np.array([0.0, -5.49, -3.39]),
               'c_wn': np.array([2670.2425, 927.92374, 831.28619]),
               'a': np.array([1.6863857, 0.39419031, 0.26364620]),
               'b': np.array([1.0 / 1.0025955,
                              1.0 / 1.0013299,
                              1.0 / 1.0009546]),
               'b0': np.array([0.0, 5.70, 3.58]),
               'b1': np.array([1 - 0.0, 1 - 0.11187, 1 - 0.05991]),
               'b2': np.array([0.0, 0.00054668, 0.00024985])}}


class Calibrator(object):

    def __init__(self, spacecraft):
        self.ah = None
        self.al = None
        self.bh = None
        self.bl = None
        self.ch = None
        self.cl = None
        self.c_s = None
        self.c_dark = None
        self.l_date = None
        self.d = None
        self.n_s = None
        self.c_wn = None
        self.a = None
        self.b = None
        self.b0 = None
        self.b1 = None
        self.b2 = None
        self.__dict__.update(coeffs[spacecraft])


def calibrate_solar(counts, chan, year, jday, spacecraft, corr=1):
    """Do the solar calibration and return reflectance (between 0 and 100)."""
    cal = Calibrator(spacecraft)

    t = (year + jday / 365.0) - cal.l_date
    stl = (cal.al[chan] * (100.0 + cal.bl[chan] * t
                           + cal.cl[chan] * t * t)) / 100.0
    sth = (cal.ah[chan] * (100.0 + cal.bh[chan] * t
                           + cal.ch[chan] * t * t)) / 100.0
    if cal.c_s is not None:
        return np.where(counts <= cal.c_s[chan],
                        (counts - cal.c_dark[chan]) * stl * corr,
                        ((cal.c_s[chan] - cal.c_dark[chan]) * stl
                         + (counts - cal.c_s[chan]) * sth) * corr)
    else:
        return (counts - cal.c_dark[chan]) * stl * corr


def calibrate_thermal(counts, prt, ict, space, line_numbers, channel, spacecraft):
    """Do the thermal calibration and return brightness temperatures (K)."""
    cal = Calibrator(spacecraft)

    chan = channel - 3

    lines, columns = counts.shape[:2]

    offset = 0

    for i, prt_val in enumerate(prt):
        if prt_val < 50:
            offset = i
            break

    iprt = (line_numbers - line_numbers[0] + 5 - offset) % 5

    ifix = np.where(np.logical_and(iprt == 1, prt < 50))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 1, prt > 50))
        prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])

    ifix = np.where(np.logical_and(iprt == 2, prt < 50))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 2, prt > 50))
        prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])

    ifix = np.where(np.logical_and(iprt == 3, prt < 50))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 3, prt > 50))
        prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])

    ifix = np.where(np.logical_and(iprt == 4, prt < 50))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 4, prt > 50))
        prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])

    tprt = (cal.d[iprt, 0] + prt *
            (cal.d[iprt, 1] + prt *
             (cal.d[iprt, 2] + prt *
              (cal.d[iprt, 3] + prt *
               (cal.d[iprt, 4])))))

    zeros = iprt == 0
    nonzeros = np.logical_not(zeros)

    tprt[zeros] = np.interp((zeros).nonzero()[0],
                            (nonzeros).nonzero()[0],
                            tprt[nonzeros])

    if channel == 3:
        zeros = ict < 100
        nonzeros = np.logical_not(zeros)

        ict[zeros] = np.interp((zeros).nonzero()[0],
                               (nonzeros).nonzero()[0],
                               ict[nonzeros])
        zeros = space < 100
        nonzeros = np.logical_not(zeros)

        space[zeros] = np.interp((zeros).nonzero()[0],
                                 (nonzeros).nonzero()[0],
                                 space[nonzeros])

    # convolving and smoothing PRT, ICT and SPACE values
    if lines > 51:
        wlength = 51
    else:
        wlength = 3

    weighting_function = np.ones(wlength, dtype=float) / wlength
    tprt_convolved = np.convolve(tprt, weighting_function, 'same')
    ict_convolved = np.convolve(ict, weighting_function, 'same')
    space_convolved = np.convolve(space, weighting_function, 'same')

    # take care of the beginning and end
    tprt_convolved[0:(wlength - 1) // 2] = tprt_convolved[(wlength - 1) // 2]
    ict_convolved[0:(wlength - 1) // 2] = ict_convolved[(wlength - 1) // 2]
    space_convolved[0:(wlength - 1) // 2] = space_convolved[(wlength - 1) // 2]
    tprt_convolved[-(wlength - 1) // 2:] = tprt_convolved[-((wlength + 1) // 2)]
    ict_convolved[-(wlength - 1) // 2:] = ict_convolved[-((wlength + 1) // 2)]
    space_convolved[-(wlength - 1) // 2:] = \
        space_convolved[-((wlength + 1) // 2)]

    new_tprt = np.transpose(np.tile(tprt_convolved, (columns, 1)))
    new_ict = np.transpose(np.tile(ict_convolved, (columns, 1)))
    new_space = np.transpose(np.tile(space_convolved, (columns, 1)))

    # calibrating thermal channel

    tBB = new_tprt
    tsBB = cal.a[chan] + cal.b[chan] * tBB
    nBB_num = (1.1910427 * 0.000010) * cal.c_wn[chan] ** 3
    nBB = nBB_num / (np.exp((1.4387752 * cal.c_wn[chan]) / tsBB) - 1.0)

    Nlin = (cal.n_s[chan] +
            (((nBB - cal.n_s[chan])
              * (new_space - counts.astype(float)))
             / (new_space - new_ict)))
    Ncor = cal.b0[chan] + Nlin * (cal.b1[chan] + cal.b2[chan] * Nlin)
    Ne = Ncor
    tsE = ((1.4387752 * cal.c_wn[chan])
           / np.log(1.0 + nBB_num / Ne))
    bt = (tsE - cal.a[chan]) / cal.b[chan]

    if chan == 0:
        bt = np.where((counts - new_space) >= 0, 0.0, bt)

    return bt
