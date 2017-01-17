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
from ._filter import _mean_filter


def mean_filter(data, fill_value, box_size):
    """
    Filter the given 2D array 'data' by averaging the valid elements within a
    box of size (boxsize x boxsize) around each pixel. Fill values are not
    taken into account.

    @param data: 2D array to be filtered
    @param box_size: Specifies the boxsize. Must be odd.
    @param fill_value: Value indicating invalid/missing data

    @return: The filtered array.
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
    """
    Filter the given 2D array 'data' by computing the standard deviation of the
    valid elements within a box of size (box_size x box_size) around each pixel.

    Masked values are not taken into account. Since

        std = sqrt( mean(data^2) - mean(data)^2 )

    we can use mean_filter() to compute the standard deviation.

    @param data: 2D array to be filtered
    @param fill_value: Value indicating invalid/missing data
    @param box_size: Specifies the boxsize. Must be odd.
    @return: The filtered array
    """
    mean_squared = np.square(mean_filter(data, box_size=box_size,
                                         fill_value=fill_value))
    squared_mean = mean_filter(np.square(data), box_size=box_size,
                               fill_value=fill_value)
    return np.ma.sqrt(squared_mean - mean_squared)


def get_tsm_idx(tar1, tar2, tar4, tar5): 
    """
    Return index of pixels where TSM issue occurs.
    """

    # absolute difference because ch1 is very similar to ch2
    abs_d12 = abs(tar1 - tar2)
    # relative difference because ch4 and ch5 differ
    rel_d45 = 100.0*(tar4 - tar5)/tar5

    # standard deviation of abs_d12 and rel_d45
    box_size = 3
    fill_value = -9999.0
    std_d12 = std_filter(abs_d12, box_size, fill_value)
    std_d45 = std_filter(rel_d45, box_size, fill_value)

    # using ch1, ch2, ch4, ch5 in combination
    # all channels seems to be affected throughout the whole orbit,
    # independent of VIS and NIR or day and night
    idx = np.where( (std_d12 > 0.02) & (std_d45 > 2.00) )

    return idx


def flag_pixels(channel1, channel2, channel3b,
                channel4, channel5, channel3a, fillv):
    """
    Set bad pixels due to temporary scan motor issue to fill value.
    Scale reflectances ranging from 0 to 1.5
    Scale brightness temperature ranging from 170 to 350.
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
        array[:,:] = np.ma.filled(array, fillv)

    # ------------------------------------------------------------------------
    # (3) Re-scaling measurments
    # ------------------------------------------------------------------------
    # re-scaling reflectance obs
    for array in [ch1, ch2, ch3a]:
        if np.ma.count(array[array!=fillv]) > 0: 
            array[array!=fillv] = (array[array!=fillv] - ref_offs) / ref_gain
    # re-scaling brightness temperature obs
    for array in [ch3b, ch4, ch5]:
        if np.ma.count(array[array!=fillv]) > 0: 
            array[array!=fillv] = (array[array!=fillv] - bt_offs) / bt_gain

    return ch1, ch2, ch3b, ch4, ch5, ch3a

