#!/usr/bin/env python

# Copyright (c) 2014-2015, 2019 Pytroll Developers

# Author(s):

#   Jonathan Mittaz <j.mittaz@reading.ac.uk>

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

"""Uncertainty information based on the noaa.py calibration for IR channels
"""
from __future__ import division

import numpy as np
import xarray as xr

from pygac.uncertainty.ir import ir_uncertainty
from pygac.uncertainty.vis import vis_uncertainty


def uncertainty(ds, mask):
    """Get and merge uncertainties from the visible and IR channels."""

    irdata = ir_uncertainty(ds,mask)
    visdata = vis_uncertainty(ds,mask)

    #
    # Get required output size (3a/3b present)
    #
    nb_refl_channels = 2
    if "3a" in ds["channels"]:
        nb_refl_channels = 3

    #
    # Merge IR/Vis uncertainties
    #
    random = np.empty_like(ds["channels"], dtype=np.float32)
    systematic = np.empty_like(ds["channels"], dtype=np.float32)

    random[:, :, 0:nb_refl_channels] = visdata["random"].values[:, :, 0:nb_refl_channels]
    random[:, :, nb_refl_channels:] = irdata["random"].values[:,:,:]
    systematic[:, :, 0:nb_refl_channels] = visdata["systematic"].values[:, :, 0:nb_refl_channels]
    systematic[:, :, nb_refl_channels:] = irdata["systematic"].values[:, :, :]

    #
    # Make xarray data arrays
    #
    random_da = xr.DataArray(random, dims=["scan_line_index","columns","channel_names"],
                             attrs={"long_name":"Random uncertainties", "units":"Albedo/K"})
    sys_da = xr.DataArray(systematic, dims=["scan_line_index","columns","channel_names"],
                          attrs={"long_name":"Systematic uncertainties", "units":"Albedo/K"})

    uratio_da = xr.DataArray(irdata["chan_covar_ratio"].values,
                             dims=["scan_line_index","columns","ir_channel_names"],
                             attrs={"long_name":"Channel-to-channel covariance  ratio"})

    #
    # Now merge flags
    #
    uflags = np.zeros_like(ds["latitude"], dtype=np.int8)
    for i in range(uflags.shape[0]):
        # IR flags are per scanline
        uflags[i,:] = irdata["uncert_flags"].values[i]
        # Vis flags are at pixel level - set 4th bit for solar contamination
        gd = (visdata["solar_fov_contam"].values[i,:] == 1)
        uflags[i,gd] = (uflags[i,gd]|8)

    uflags_da = xr.DataArray(uflags, dims=["scan_line_index","columns"],
                             attrs={"long_name": "Uncertainty flags",
                                    "flag_masks": "1b, 2b, 4b, 8b",
                                    "flag_meanings": ("bad_space_view "
                                                      "solar_contamination_of_gain "
                                                      "no_IR_systematic_uncertainty "
                                                      "solar_contamination_of_FOV")})

    uncertainties = xr.Dataset(dict(random=random_da,systematic=sys_da,
                                    chan_covar_ratio=uratio_da,
                                    uncert_flags=uflags_da))

    return uncertainties
