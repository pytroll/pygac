#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author(s):

#   Stephan Finkensieper <stephan.finkensieper@dwd.de>

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

"""Test TSM module."""


import unittest
import pygac.correct_tsm_issue as tsm

import numpy as np
import numpy.testing


class TSMTest(unittest.TestCase):
    """Test TSM module."""

    def test_get_tsm_idx(self):
        """Test identification of TSM affected pixels."""
        cols = np.arange(409)
        rows = np.arange(3000)
        mcols, mrows = np.meshgrid(cols, rows)

        # Create dummy reflectances in [0, 1] and BTs [170, 350]
        refl = np.sin(mrows / float(rows.size)) * \
            np.cos(mcols / float(cols.size))
        bt = np.cos(mrows / float(rows.size)) * \
            np.sin(mcols / float(cols.size))
        bt = 170 + (350 - 170) * bt

        ch1 = refl.copy()
        ch2 = refl.copy()
        ch4 = bt.copy()
        ch5 = bt.copy()

        # Add rectangle with noise to channel 1 and 4
        noise_cols, noise_rows = np.meshgrid(np.arange(200, 300),
                                             np.arange(1000, 2000))
        for ch in [ch1, ch4]:
            ch[noise_rows, noise_cols] += 1000 + 1000*np.random.rand(
                noise_rows.shape[0], noise_rows.shape[1])

        # We expect the filter to also detect the edges of the noisy
        # rectangle: At the edges there is a transition from zero to
        # nonzero channel diff which leads to an increased 3x3 standard
        # deviation.
        noise_cols_exp, noise_rows_exp = np.meshgrid(np.arange(199, 301),
                                                     np.arange(999, 2001))
        idx = tsm.get_tsm_idx(ch1, ch2, ch4, ch5)
        numpy.testing.assert_array_equal(idx[0], noise_rows_exp.ravel())
        numpy.testing.assert_array_equal(idx[1], noise_cols_exp.ravel())

    def test_mean_filter(self):
        """Test gridbox mean filter."""
        # Define test data
        data = np.array(
            [[1, 2, 2, 1],
             [2, np.nan, np.nan, 1],
             [1, 1, 2, 2],
             [2, 2, 1, np.nan]]
        )

        # Define reference
        filtered_ref = np.array(
            [[5/3., 7/4., 6/4., 4/3.],
             [7/5., 11/7., 11/7., 8/5.],
             [8/5., 11/7., 9/6., 6/4.],
             [6/4., 9/6., 8/5., 5/3.]]
        )

        # Apply mean and compare results against reference
        filtered = tsm.mean_filter(data=data, box_size=3, fill_value=-999)
        numpy.testing.assert_array_equal(filtered, filtered_ref)


def suite():
    """The suite for test_corrections"""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TSMTest))

    return mysuite


if __name__ == '__main__':
    unittest.main()
