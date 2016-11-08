#
# Copyright (c) Nov 2016, DWD
#
# Author(s):
#
#   Stephan Finkensieper <stephan.finkensieper@dwd.dwd>
#   Cornelia Schlundt <cornelia.schlundt@dwd.de>
#
# Purpose
#   NOAA-14, 15 and 16 suffer from temporary scan motor issue
#   i.e. parts of an orbit contain corrupt data,
#        which are flagged with fill_values
#

import os
import numpy as np
import weave


def gridbox_mean(data, fill_value, box_size):
    """
    For each element in C{data}, calculate the mean value of all valid elements
    in the C{box_size} x C{box_size} box around it.

    Masked values are not taken into account.

    @param data: 2D masked array
    @param fill_value: Value to be used to fill masked elements
    @param box_size: Box size
    @return: Gridbox mean
    @rtype: numpy.ma.core.MaskedArray
    """
    if not box_size % 2 == 1:
        raise ValueError('Box size must be odd.')

    # Fill masked elements
    fdata = data.astype('f8').filled(fill_value)

    # Allocate filtered image
    filtered = np.zeros(data.shape, dtype='f8')
    nrows, ncols = data.shape

    c_code = """
    int row, col, rowbox, colbox, nbox;
    double fill_value_d = (double) fill_value;
    int radius = (box_size-1)/2;


    for(row=0; row<nrows; row++)
    {
        for(col=0; col<ncols; col++)
        {
            filtered(row,col) = 0;
            nbox = 0;
            for(rowbox=row-radius; rowbox<=row+radius; rowbox++)
            {
                for(colbox=col-radius; colbox<=col+radius; colbox++)
                {
                    if(rowbox >= 0 && rowbox < nrows && colbox >= 0 and colbox < ncols)
                    {
                        if(fdata(rowbox,colbox) != fill_value_d)
                        {
                            filtered(row,col) += fdata(rowbox,colbox);
                            nbox += 1;
                        }
                    }
                }
            }
            if(nbox > 0)
            {
                filtered(row,col) /= (double) nbox;
            }
        }
    }
    return_val=0;
    """

    # Execute inline C code
    err = weave.inline(
        c_code,
        ['fdata', 'filtered', 'fill_value', 'box_size', 'nrows', 'ncols'],
        type_converters=weave.converters.blitz,
        compiler="gcc"
    )
    if err != 0:
        raise RuntimeError('Blitz failed with returncode {0}'.format(err))

    # Re-mask fill values
    return np.ma.masked_equal(filtered, fill_value)


def gridbox_std(data, box_size, fill_value):
    """
    For each element in C{data}, calculate the standard deviation of all valid
    elements in the C{box_size} x C{box_size} box around it.

    Masked values are not taken into account. Since

        std = sqrt( mean(data^2) - mean(data)^2 )

    we can use L{gridbox_mean} to compute the standard deviation.

    @param data: 2D masked array
    @param fill_value: Value to be used to fill masked elements
    @param box_size: Box size
    @return: Gridbox standard deviation
    @rtype: numpy.ma.core.MaskedArray
    """
    mean_squared = np.square(gridbox_mean(data, box_size=box_size, fill_value=fill_value))
    squared_mean = gridbox_mean(np.square(data), box_size=box_size, fill_value=fill_value)
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
    std_d12 = gridbox_std(abs_d12, box_size, fill_value)
    std_d45 = gridbox_std(rel_d45, box_size, fill_value)

    # using ch1, ch2, ch4, ch5 in combination
    # all channels seems to be affected throughout the whole orbit,
    # independent of VIS and NIR or day and night
    idx = np.where( (std_d12 > 0.02) & (std_d45 > 2.00) )

    return idx


def flag_pixels(c1, c2, c3, c4, c5, c6, fillv):
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
    ch1 = ref_gain * np.ma.masked_equal(c1, fillv) + ref_offs
    ch2 = ref_gain * np.ma.masked_equal(c2, fillv) + ref_offs
    ch6 = ref_gain * np.ma.masked_equal(c6, fillv) + ref_offs
    # original bt data
    ch3 = bt_gain * np.ma.masked_equal(c3, fillv) + bt_offs
    ch4 = bt_gain * np.ma.masked_equal(c4, fillv) + bt_offs
    ch5 = bt_gain * np.ma.masked_equal(c5, fillv) + bt_offs
    
    # ------------------------------------------------------------------------
    # (2) TSM Correction
    # ------------------------------------------------------------------------
    # find indices of tsm issue affected pixels
    idx = get_tsm_idx(ch1, ch2, ch4, ch5)
    # apply correction index using fill_value and fill masked elements
    for array in [ch1, ch2, ch3, ch4, ch5, ch6]:
        if isinstance(array.mask, np.bool_): 
            array.mask = np.zeros(array.shape, dtype='bool') 
        array.mask[idx] = True
        array[:,:] = np.ma.filled(array, fillv)

    # ------------------------------------------------------------------------
    # (3) Re-scaling measurments
    # ------------------------------------------------------------------------
    # re-scaling reflectance obs
    for array in [ch1, ch2, ch6]:
        if np.ma.count(array[array!=fillv]) > 0: 
            array[array!=fillv] = (array[array!=fillv] - ref_offs) / ref_gain
    # re-scaling brightness temperature obs
    for array in [ch3, ch4, ch5]:
        if np.ma.count(array[array!=fillv]) > 0: 
            array[array!=fillv] = (array[array!=fillv] - bt_offs) / bt_gain

    return ch1, ch2, ch3, ch4, ch5, ch6

