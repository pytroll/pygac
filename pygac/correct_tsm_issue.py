#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author(s):

#   Stephan Finkensieper <stephan.finkensieper@dwd.dwd>
#   Cornelia Schlundt <cornelia.schlundt@dwd.de>

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

"""NOAA-14, 15 and 16 suffer from temporary scan motor issue i.e. parts of
an orbit contain corrupt data. This module identifies affected pixels and flags
them with fill values."""

import numpy as np
import datetime
from ._filter import _mean_filter

TSM_AFFECTED_INTERVALS_POD = {
    3: [(datetime.datetime(2001, 10, 19, 4, 50), datetime.datetime(2001, 10, 19, 13, 38)),
        (datetime.datetime(2001, 10, 19, 16, 58), datetime.datetime(2001, 10, 19, 18, 7)),
        (datetime.datetime(2001, 10, 19, 23, 55), datetime.datetime(2001, 10, 27, 3, 19)),
        (datetime.datetime(2001, 10, 27, 4, 53), datetime.datetime(2001, 10, 31, 21, 7)),
        (datetime.datetime(2001, 10, 31, 20, 58), datetime.datetime(2001, 11, 13, 12, 46)),
        (datetime.datetime(2001, 11, 13, 13, 25), datetime.datetime(2001, 11, 21, 23, 40)),
        (datetime.datetime(2001, 11, 22, 1, 14), datetime.datetime(2001, 11, 25, 6, 1)),
        (datetime.datetime(2001, 11, 25, 11, 3), datetime.datetime(2001, 11, 25, 13, 8)),
        (datetime.datetime(2001, 11, 26, 10, 51), datetime.datetime(2001, 11, 26, 12, 55)),
        (datetime.datetime(2001, 11, 26, 15, 49), datetime.datetime(2001, 11, 26, 19, 26)),
        (datetime.datetime(2001, 11, 26, 20, 36), datetime.datetime(2001, 11, 26, 22, 31)),
        (datetime.datetime(2001, 11, 27, 0, 5), datetime.datetime(2001, 12, 2, 23, 7)),
        (datetime.datetime(2001, 12, 3, 0, 41), datetime.datetime(2001, 12, 15, 0, 2)),
        (datetime.datetime(2001, 12, 15, 1, 30), datetime.datetime(2001, 12, 15, 3, 26)),
        (datetime.datetime(2001, 12, 15, 5, 2), datetime.datetime(2001, 12, 15, 10, 36)),
        (datetime.datetime(2001, 12, 15, 12, 2), datetime.datetime(2001, 12, 15, 13, 53)),
        (datetime.datetime(2001, 12, 15, 15, 21), datetime.datetime(2001, 12, 15, 18, 56)),
        (datetime.datetime(2001, 12, 15, 23, 37), datetime.datetime(2001, 12, 16, 1, 33)),
        (datetime.datetime(2002, 1, 1, 10, 14), datetime.datetime(2002, 1, 10, 1, 29)),
        (datetime.datetime(2002, 7, 26, 14, 40), datetime.datetime(2002, 7, 27, 13, 6)),
        (datetime.datetime(2002, 7, 29, 5, 26), datetime.datetime(2002, 8, 8, 0, 2)),
        (datetime.datetime(2002, 8, 8, 1, 32), datetime.datetime(2002, 8, 12, 22, 51)),
        (datetime.datetime(2002, 8, 13, 0, 25), datetime.datetime(2002, 9, 30, 4, 8))
        ]  # NOAA-14
}
TSM_AFFECTED_INTERVALS_KLM = {
    4: [(datetime.datetime(2001, 1, 8, 23, 55), datetime.datetime(2001, 1, 9, 23, 47)),
        (datetime.datetime(2001, 1, 10, 1, 14), datetime.datetime(2001, 1, 11, 1, 6)),
        (datetime.datetime(2001, 1, 23, 18, 32), datetime.datetime(2001, 1, 26, 0, 25)),
        (datetime.datetime(2001, 1, 30, 10, 59), datetime.datetime(2001, 1, 30, 17, 58)),
        (datetime.datetime(2001, 2, 2, 2, 42), datetime.datetime(2001, 2, 8, 11, 10)),
        (datetime.datetime(2001, 2, 11, 6, 3), datetime.datetime(2001, 2, 16, 21, 32)),
        (datetime.datetime(2001, 4, 24, 5, 25), datetime.datetime(2001, 4, 24, 9, 18)),
        (datetime.datetime(2001, 4, 24, 10, 53), datetime.datetime(2001, 4, 25, 8, 55)),
        (datetime.datetime(2001, 4, 26, 2, 59), datetime.datetime(2001, 4, 26, 8, 25)),
        (datetime.datetime(2001, 4, 26, 18, 22), datetime.datetime(2001, 4, 26, 20, 18)),
        (datetime.datetime(2001, 4, 26, 23, 5), datetime.datetime(2001, 4, 27, 8, 2)),
        (datetime.datetime(2001, 4, 28, 22, 23), datetime.datetime(2001, 4, 29, 0, 14)),
        (datetime.datetime(2001, 4, 29, 15, 33), datetime.datetime(2001, 4, 29, 17, 28)),
        (datetime.datetime(2001, 4, 29, 20, 33), datetime.datetime(2001, 4, 30, 8, 35)),
        (datetime.datetime(2001, 4, 30, 13, 35), datetime.datetime(2001, 4, 30, 18, 46)),
        (datetime.datetime(2001, 4, 30, 20, 11), datetime.datetime(2001, 5, 1, 8, 11)),
        (datetime.datetime(2001, 5, 2, 2, 23), datetime.datetime(2001, 5, 2, 5, 58)),
        (datetime.datetime(2001, 5, 2, 9, 22), datetime.datetime(2001, 5, 2, 22, 45)),
        (datetime.datetime(2001, 5, 3, 0, 8), datetime.datetime(2001, 5, 3, 9, 14)),
        (datetime.datetime(2001, 5, 8, 10, 34), datetime.datetime(2001, 5, 8, 14, 8)),
        (datetime.datetime(2001, 5, 8, 15, 30), datetime.datetime(2001, 5, 9, 15, 34)),
        (datetime.datetime(2001, 5, 10, 0, 50), datetime.datetime(2001, 5, 10, 4, 35)),
        (datetime.datetime(2001, 5, 10, 13, 7), datetime.datetime(2001, 5, 10, 21, 48)),
        (datetime.datetime(2001, 5, 10, 22, 47), datetime.datetime(2001, 5, 11, 0, 42)),
        (datetime.datetime(2001, 5, 11, 3, 57), datetime.datetime(2001, 5, 11, 7, 43)),
        (datetime.datetime(2001, 5, 11, 11, 5), datetime.datetime(2001, 5, 11, 17, 55)),
        (datetime.datetime(2001, 5, 11, 19, 21), datetime.datetime(2001, 5, 12, 7, 19)),
        (datetime.datetime(2001, 5, 12, 8, 56), datetime.datetime(2001, 5, 14, 6, 26)),
        (datetime.datetime(2001, 5, 14, 21, 32), datetime.datetime(2001, 5, 14, 23, 11)),
        (datetime.datetime(2001, 5, 15, 4, 6), datetime.datetime(2001, 5, 15, 7, 52)),
        (datetime.datetime(2001, 5, 15, 17, 59), datetime.datetime(2001, 5, 15, 19, 44)),
        (datetime.datetime(2001, 5, 15, 21, 19), datetime.datetime(2001, 5, 16, 7, 29)),
        (datetime.datetime(2001, 5, 16, 17, 26), datetime.datetime(2001, 5, 16, 19, 22)),
        (datetime.datetime(2001, 5, 16, 23, 50), datetime.datetime(2001, 5, 17, 8, 55)),
        (datetime.datetime(2001, 5, 17, 15, 26), datetime.datetime(2001, 5, 17, 17, 18)),
        (datetime.datetime(2001, 5, 17, 20, 25), datetime.datetime(2001, 5, 18, 6, 36)),
        (datetime.datetime(2001, 6, 5, 16, 34), datetime.datetime(2001, 6, 5, 21, 46)),
        (datetime.datetime(2001, 6, 6, 0, 35), datetime.datetime(2001, 6, 6, 7, 51)),
        (datetime.datetime(2001, 8, 4, 5, 2), datetime.datetime(2001, 8, 4, 8, 54)),
        (datetime.datetime(2001, 8, 6, 21, 17), datetime.datetime(2001, 8, 7, 7, 37)),
        (datetime.datetime(2001, 8, 7, 19, 16), datetime.datetime(2001, 8, 8, 12, 31)),
        (datetime.datetime(2001, 8, 8, 15, 33), datetime.datetime(2001, 8, 8, 22, 13)),
        (datetime.datetime(2001, 8, 23, 0, 58), datetime.datetime(2001, 8, 23, 15, 18)),
        (datetime.datetime(2001, 8, 23, 21, 30), datetime.datetime(2001, 8, 24, 7, 50)),
        (datetime.datetime(2001, 8, 24, 9, 24), datetime.datetime(2001, 8, 25, 0, 25)),
        (datetime.datetime(2001, 8, 25, 3, 50), datetime.datetime(2001, 8, 25, 11, 5)),
        (datetime.datetime(2001, 8, 25, 12, 28), datetime.datetime(2001, 8, 25, 17, 39)),
        (datetime.datetime(2001, 8, 25, 20, 54), datetime.datetime(2001, 8, 26, 8, 52)),
        (datetime.datetime(2001, 8, 26, 18, 41), datetime.datetime(2001, 8, 27, 4, 51)),
        (datetime.datetime(2001, 8, 27, 23, 2), datetime.datetime(2001, 8, 28, 0, 58)),
        (datetime.datetime(2001, 8, 28, 4, 14), datetime.datetime(2001, 8, 28, 7, 59)),
        (datetime.datetime(2001, 8, 28, 14, 49), datetime.datetime(2001, 8, 28, 16, 31)),
        (datetime.datetime(2001, 8, 28, 17, 58), datetime.datetime(2001, 8, 29, 17, 48)),
        (datetime.datetime(2001, 8, 29, 22, 19), datetime.datetime(2001, 8, 30, 0, 11)),
        (datetime.datetime(2001, 8, 30, 5, 7), datetime.datetime(2001, 8, 31, 13, 46)),
        (datetime.datetime(2001, 8, 31, 18, 37), datetime.datetime(2001, 8, 31, 22, 11)),
        (datetime.datetime(2001, 8, 31, 23, 10), datetime.datetime(2001, 9, 1, 8, 8)),
        (datetime.datetime(2001, 9, 1, 11, 29), datetime.datetime(2001, 9, 1, 13, 23)),
        (datetime.datetime(2001, 9, 1, 14, 58), datetime.datetime(2001, 9, 1, 16, 39)),
        (datetime.datetime(2001, 9, 1, 18, 13), datetime.datetime(2001, 9, 1, 23, 3)),
        (datetime.datetime(2001, 9, 2, 0, 28), datetime.datetime(2001, 9, 2, 7, 44)),
        (datetime.datetime(2001, 9, 2, 16, 5), datetime.datetime(2001, 9, 2, 18, 8)),
        (datetime.datetime(2001, 9, 2, 22, 26), datetime.datetime(2001, 9, 3, 9, 10)),
        (datetime.datetime(2001, 9, 3, 12, 22), datetime.datetime(2001, 9, 4, 5, 8)),
        (datetime.datetime(2001, 9, 7, 2, 4), datetime.datetime(2001, 9, 7, 9, 19)),
        (datetime.datetime(2001, 9, 7, 10, 54), datetime.datetime(2001, 9, 7, 19, 23)),
        (datetime.datetime(2001, 9, 7, 23, 51), datetime.datetime(2001, 9, 8, 20, 40)),
        (datetime.datetime(2001, 9, 8, 23, 29), datetime.datetime(2001, 9, 10, 22, 58)),
        (datetime.datetime(2001, 9, 11, 0, 25), datetime.datetime(2001, 9, 11, 4, 16)),
        (datetime.datetime(2001, 9, 11, 5, 35), datetime.datetime(2001, 9, 11, 12, 55)),
        (datetime.datetime(2001, 10, 25, 17, 39), datetime.datetime(2001, 10, 26, 0, 8)),
        (datetime.datetime(2002, 3, 17, 14, 29), datetime.datetime(2002, 3, 18, 2, 4)),
        (datetime.datetime(2002, 3, 18, 5, 11), datetime.datetime(2002, 3, 18, 7, 15)),
        (datetime.datetime(2002, 3, 18, 18, 53), datetime.datetime(2002, 3, 18, 20, 48))
        ],   # NOAA-15
    2: [(datetime.datetime(2004, 1, 14, 14, 13), datetime.datetime(2004, 1, 14, 22, 46)),
        (datetime.datetime(2004, 1, 15, 15, 45), datetime.datetime(2004, 1, 15, 19, 3)),
        (datetime.datetime(2004, 1, 15, 22, 28), datetime.datetime(2004, 1, 16, 2, 6)),
        (datetime.datetime(2004, 1, 16, 5, 29), datetime.datetime(2004, 1, 16, 22, 23)),
        (datetime.datetime(2004, 1, 16, 23, 58), datetime.datetime(2004, 1, 18, 0, 1)),
        (datetime.datetime(2004, 1, 18, 1, 23), datetime.datetime(2004, 1, 18, 8, 45)),
        (datetime.datetime(2004, 1, 18, 10, 6), datetime.datetime(2004, 1, 18, 23, 42)),
        (datetime.datetime(2004, 1, 19, 1, 11), datetime.datetime(2004, 1, 19, 15, 12)),
        (datetime.datetime(2004, 1, 19, 19, 44), datetime.datetime(2004, 1, 19, 23, 30)),
        (datetime.datetime(2004, 1, 20, 0, 59), datetime.datetime(2004, 1, 20, 4, 52)),
        (datetime.datetime(2004, 1, 20, 9, 56), datetime.datetime(2004, 1, 20, 15, 1)),
        (datetime.datetime(2004, 1, 20, 16, 37), datetime.datetime(2004, 1, 21, 4, 40)),
        (datetime.datetime(2004, 1, 21, 11, 13), datetime.datetime(2004, 1, 21, 16, 39)),
        (datetime.datetime(2004, 1, 22, 0, 35), datetime.datetime(2004, 1, 22, 6, 18)),
        (datetime.datetime(2004, 1, 22, 12, 42), datetime.datetime(2004, 1, 23, 2, 27)),
        (datetime.datetime(2004, 1, 23, 10, 50), datetime.datetime(2004, 1, 23, 22, 45)),
        (datetime.datetime(2004, 1, 24, 0, 12), datetime.datetime(2004, 1, 24, 5, 54)),
        (datetime.datetime(2004, 3, 16, 22, 38), datetime.datetime(2004, 3, 17, 22, 41)),
        (datetime.datetime(2004, 3, 18, 0, 1), datetime.datetime(2004, 3, 18, 5, 44)),
        (datetime.datetime(2004, 3, 20, 18, 13), datetime.datetime(2004, 3, 20, 23, 42)),
        (datetime.datetime(2004, 3, 21, 2, 59), datetime.datetime(2004, 3, 21, 5, 3)),
        (datetime.datetime(2004, 3, 21, 10, 9), datetime.datetime(2004, 3, 21, 23, 29)),
        (datetime.datetime(2004, 3, 22, 0, 58), datetime.datetime(2004, 3, 22, 6, 41)),
        (datetime.datetime(2004, 3, 22, 9, 55), datetime.datetime(2004, 3, 22, 15, 0)),
        (datetime.datetime(2004, 3, 22, 16, 26), datetime.datetime(2004, 3, 22, 23, 18)),
        (datetime.datetime(2004, 3, 23, 0, 46), datetime.datetime(2004, 3, 23, 4, 39)),
        (datetime.datetime(2004, 3, 23, 6, 14), datetime.datetime(2004, 3, 23, 11, 27)),
        (datetime.datetime(2004, 3, 23, 12, 53), datetime.datetime(2004, 3, 23, 16, 38)),
        (datetime.datetime(2004, 3, 23, 19, 20), datetime.datetime(2004, 3, 24, 6, 17)),
        (datetime.datetime(2004, 3, 24, 7, 43), datetime.datetime(2004, 3, 25, 6, 5)),
        (datetime.datetime(2004, 3, 25, 20, 46), datetime.datetime(2004, 3, 25, 22, 50)),
        (datetime.datetime(2004, 3, 26, 0, 11), datetime.datetime(2004, 3, 26, 5, 54)),
        (datetime.datetime(2004, 3, 26, 10, 39), datetime.datetime(2004, 3, 26, 14, 15)),
        (datetime.datetime(2004, 3, 28, 8, 40), datetime.datetime(2004, 4, 1, 23, 5)),
        (datetime.datetime(2004, 4, 2, 0, 32), datetime.datetime(2004, 4, 2, 22, 53)),
        (datetime.datetime(2004, 4, 3, 0, 21), datetime.datetime(2004, 4, 4, 22, 31)),
        (datetime.datetime(2004, 4, 5, 0, 6), datetime.datetime(2004, 4, 5, 22, 19)),
        (datetime.datetime(2004, 4, 5, 23, 55), datetime.datetime(2004, 4, 6, 5, 23)),
        (datetime.datetime(2004, 4, 6, 18, 21), datetime.datetime(2004, 4, 6, 23, 49)),
        (datetime.datetime(2004, 4, 7, 13, 24), datetime.datetime(2004, 4, 7, 15, 20)),
        (datetime.datetime(2004, 4, 7, 16, 45), datetime.datetime(2004, 4, 7, 23, 38)),
        (datetime.datetime(2004, 4, 8, 1, 6), datetime.datetime(2004, 4, 9, 10, 7)),
        (datetime.datetime(2004, 4, 21, 2, 1), datetime.datetime(2004, 4, 21, 17, 41)),
        (datetime.datetime(2004, 4, 21, 18, 51), datetime.datetime(2004, 4, 21, 22, 38)),
        (datetime.datetime(2004, 4, 22, 0, 13), datetime.datetime(2004, 4, 22, 5, 43)),
        (datetime.datetime(2004, 4, 24, 21, 49), datetime.datetime(2004, 4, 24, 23, 47)),
        (datetime.datetime(2004, 4, 25, 1, 14), datetime.datetime(2004, 4, 25, 20, 3)),
        (datetime.datetime(2004, 4, 26, 2, 51), datetime.datetime(2004, 4, 26, 11, 43)),
        (datetime.datetime(2004, 4, 26, 21, 27), datetime.datetime(2004, 4, 26, 23, 22)),
        (datetime.datetime(2004, 4, 27, 0, 50), datetime.datetime(2004, 4, 27, 4, 43)),
        (datetime.datetime(2004, 4, 27, 8, 0), datetime.datetime(2004, 4, 27, 10, 3)),
        (datetime.datetime(2004, 4, 27, 12, 57), datetime.datetime(2004, 4, 27, 21, 28)),
        (datetime.datetime(2004, 4, 28, 0, 38), datetime.datetime(2004, 4, 28, 4, 32)),
        (datetime.datetime(2004, 4, 28, 11, 5), datetime.datetime(2004, 4, 28, 16, 22)),
        (datetime.datetime(2004, 4, 29, 0, 26), datetime.datetime(2004, 4, 29, 14, 30)),
        (datetime.datetime(2004, 5, 2, 8, 44), datetime.datetime(2004, 5, 2, 23, 55)),
        (datetime.datetime(2004, 5, 3, 1, 23), datetime.datetime(2004, 5, 4, 3, 16)),
        (datetime.datetime(2004, 5, 4, 4, 50), datetime.datetime(2004, 5, 4, 8, 35)),
        (datetime.datetime(2004, 5, 4, 10, 9), datetime.datetime(2004, 5, 4, 23, 32)),
        (datetime.datetime(2004, 5, 5, 1, 0), datetime.datetime(2004, 5, 5, 3, 4)),
        (datetime.datetime(2004, 5, 5, 6, 27), datetime.datetime(2004, 5, 5, 10, 13)),
        (datetime.datetime(2004, 5, 5, 11, 27), datetime.datetime(2004, 5, 5, 13, 22)),
        (datetime.datetime(2004, 5, 8, 9, 23), datetime.datetime(2004, 5, 8, 16, 9)),
        (datetime.datetime(2004, 5, 9, 7, 22), datetime.datetime(2004, 5, 9, 10, 55)),
        (datetime.datetime(2004, 5, 10, 22, 8), datetime.datetime(2004, 5, 11, 3, 38)),
        (datetime.datetime(2004, 5, 12, 21, 44), datetime.datetime(2004, 5, 13, 3, 13)),
        (datetime.datetime(2004, 5, 13, 6, 37), datetime.datetime(2004, 5, 13, 10, 21)),
        (datetime.datetime(2004, 5, 13, 11, 35), datetime.datetime(2004, 5, 13, 15, 12)),
        (datetime.datetime(2004, 5, 13, 18, 2), datetime.datetime(2004, 5, 14, 1, 13)),
        (datetime.datetime(2004, 5, 14, 6, 25), datetime.datetime(2004, 5, 14, 11, 39)),
        (datetime.datetime(2004, 5, 14, 13, 5), datetime.datetime(2004, 5, 14, 16, 41)),
        (datetime.datetime(2004, 5, 15, 6, 14), datetime.datetime(2004, 5, 15, 16, 30)),
        (datetime.datetime(2004, 5, 15, 19, 22), datetime.datetime(2004, 5, 15, 21, 24)),
        (datetime.datetime(2004, 5, 19, 22, 9), datetime.datetime(2004, 5, 20, 22, 9)),
        (datetime.datetime(2004, 5, 21, 8, 28), datetime.datetime(2004, 5, 21, 10, 19))
        ]  # NOAA-16
}


def mean_filter(data, fill_value, box_size):
    """Filter a 2D array using an arithmetic mean kernel.

    Compute the arithmetic mean of the valid elements within a box of size
    (boxsize x boxsize) around each pixel. Masked elements are not taken into
    account.

    Args:
        data (numpy.ma.core.MaskedArray): 2D array to be filtered
        box_size (int): Specifies the boxsize. Must be odd.
        fill_value: Value to fill masked elements with. Must be outside the
            valid range of the data.

    Returns:
        numpy.ma.core.MaskedArray: The filtered array.
    """
    if not box_size % 2 == 1:
        raise ValueError('Box size must be odd.')

    # Replace masked elements with fill_value
    if isinstance(data, np.ma.core.MaskedArray):
        filled = data.filled(fill_value)
    else:
        filled = data

    # Convert data to double (this is what _mean_filter() is expecting) and
    # apply mean filter
    filtered = _mean_filter(data=filled.astype('f8'), box_size=box_size,
                            fill_value=fill_value)

    # Re-mask fill values
    return np.ma.masked_equal(filtered, fill_value)


def std_filter(data, box_size, fill_value):
    """Filter a 2D array using a standard deviation kernel.

    Compute the standard deviation of the valid elements within a box of size
    (box_size x box_size) around each pixel. Masked values are not taken into
    account. Since

        std = sqrt( mean(data^2) - mean(data)^2 )

    we can use mean_filter() to compute the standard deviation.

    Args:
        data (np.ma.core.MaskedArray): 2D array to be filtered
        box_size (int): Specifies the boxsize. Must be odd.
        fill_value: Value indicating invalid/missing data

    Returns:
        np.ma.core.MaskedArray: The filtered array
    """
    mean_squared = np.square(mean_filter(data, box_size=box_size,
                                         fill_value=fill_value))
    squared_mean = mean_filter(np.square(data), box_size=box_size,
                               fill_value=fill_value)
    return np.ma.sqrt(squared_mean - mean_squared)


def get_tsm_idx(ch1, ch2, ch4, ch5):
    """Determine indices of TSM affected pixels."""

    # absolute difference because ch1 is very similar to ch2
    abs_d12 = abs(ch1 - ch2)
    # relative difference because ch4 and ch5 differ
    rel_d45 = 100.0*(ch4 - ch5) / ch5

    # standard deviation of abs_d12 and rel_d45
    box_size = 3
    fill_value = -9999.0
    std_d12 = std_filter(abs_d12, box_size, fill_value)
    std_d45 = std_filter(rel_d45, box_size, fill_value)

    # using ch1, ch2, ch4, ch5 in combination
    # all channels seems to be affected throughout the whole orbit,
    # independent of VIS and NIR or day and night
    idx = np.where((std_d12 > 0.02) & (std_d45 > 2.00))

    return idx


def flag_pixels(channel1, channel2, channel3b,
                channel4, channel5, channel3a, fillv):
    """Set TSM affected pixels to fill value.

    Scale reflectances ranging from 0 to 1.5 and brightness temperatures
    ranging from 170 to 350.
    """
    # ------------------------------------------------------------------------
    # (1) Scaling measurements w.r.t. threshold values
    # ------------------------------------------------------------------------
    # ref between 0 and 1
    ref_gain = 0.01*0.01
    ref_offs = 0.0
    # bt between 170 and 350
    bt_gain = 0.01
    bt_offs = 273.15
    # original ref data
    ch1 = ref_gain * np.ma.masked_equal(channel1, fillv) + ref_offs
    ch2 = ref_gain * np.ma.masked_equal(channel2, fillv) + ref_offs
    ch3a = ref_gain * np.ma.masked_equal(channel3a, fillv) + ref_offs
    # original bt data
    ch3b = bt_gain * np.ma.masked_equal(channel3b, fillv) + bt_offs
    ch4 = bt_gain * np.ma.masked_equal(channel4, fillv) + bt_offs
    ch5 = bt_gain * np.ma.masked_equal(channel5, fillv) + bt_offs

    # ------------------------------------------------------------------------
    # (2) TSM Correction
    # ------------------------------------------------------------------------
    # find indices of tsm issue affected pixels
    idx = get_tsm_idx(ch1, ch2, ch4, ch5)
    # apply correction index using fill_value and fill masked elements
    for array in [ch1, ch2, ch3b, ch4, ch5, ch3a]:
        if isinstance(array.mask, np.bool_):
            array.mask = np.zeros(array.shape, dtype='bool')
        array.mask[idx] = True
        array[:, :] = np.ma.filled(array, fillv)

    # ------------------------------------------------------------------------
    # (3) Re-scaling measurments
    # ------------------------------------------------------------------------
    # re-scaling reflectance obs
    for array in [ch1, ch2, ch3a]:
        if np.ma.count(array[array != fillv]) > 0:
            array[array != fillv] = (array[array != fillv] - ref_offs) / ref_gain
    # re-scaling brightness temperature obs
    for array in [ch3b, ch4, ch5]:
        if np.ma.count(array[array != fillv]) > 0:
            array[array != fillv] = (array[array != fillv] - bt_offs) / bt_gain

    return ch1, ch2, ch3b, ch4, ch5, ch3a
