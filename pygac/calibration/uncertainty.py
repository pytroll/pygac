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
import argparse

from pygac import get_reader_class
from pygac.calibration.ir_uncertainty import ir_uncertainty, allan_deviation
from pygac.calibration.vis_uncertainty import vis_uncertainty

#
# Get amd merge uncertainties from the visible and IR channels
#
def uncertainty(ds,mask,plot=False):
    """Get combined vis/IR uncertainties"""

    irdata = ir_uncertainty(ds,mask)
    visdata = vis_uncertainty(ds,mask)

    #
    # Get required output size (3a/3b present)
    #
    fivechan = False
    if ds['channels'].values.shape[2] == 5:
        fivechan = True

    #
    # Setup for output
    #
    time = (ds["times"].values - np.datetime64("1970-01-01 00:00:00"))/\
           np.timedelta64(1,'s')
    time_da = xr.DataArray(time,dims=["times"],attrs={"long_name":"scanline time",
                                                     "units":"seconds since 1970-01-01"})
    across_da = xr.DataArray(np.arange(ds['channels'].shape[1]),dims=["across_track"])
    if fivechan:
        channels_da = xr.DataArray(np.array([1,2,3,4,5]),dims=["channels"],attrs={"long_name":"chan 1=0.6mu, 2=0.8mu, 3=3.7mu, 4=11mu, 5=12mu (NaN if not present)"})
        ir_channels_da = xr.DataArray(np.array([3,4,5]),dims=["ir_channels"],attrs={"long_name":"chan 3=3.7mu, 4=11mu, 5=12mu"})
    else:
        channels_da = xr.DataArray(np.array([1,2,3,4,5,6]),dims=["channels"],attrs={"long_name":"chan 1=0.6mu, 2=0.8mu, 3=1.6mu, 4=3.7mu, 5=11mu, 6=12mu"})
        ir_channels_da = xr.DataArray(np.array([4,5,6]),dims=["ir_channels"],attrs={"long_name":"chan 4=3.7mu, 5=11mu, 6=12mu"})

    #
    # Merge IR/Vis uncertainties
    #
    random = np.zeros(ds['channels'].values.shape,dtype=np.float32)
    systematic = np.zeros(ds['channels'].values.shape,dtype=np.float32)

    if fivechan:
        random[:,:,0:2] = visdata['random'].values[:,:,0:2]
        random[:,:,2:5] = irdata['random'].values[:,:,:]
        systematic[:,:,0:2] = visdata['systematic'].values[:,:,0:2]
        systematic[:,:,2:5] = irdata['systematic'].values[:,:,:]
    else:
        random[:,:,0:3] = visdata['random'].values[:,:,0:3]
        random[:,:,3:6] = irdata['random'].values[:,:,:]
        systematic[:,:,0:3] = visdata['systematic'].values[:,:,0:3]
        systematic[:,:,3:6] = irdata['systematic'].values[:,:,:]
                                
    #
    # Make xarray data arrays
    #
    random_da = xr.DataArray(random,dims=["times","across_track","channels"],
                             attrs={"long_name":"Random uncertainties","units":"Albedo/K"})
    sys_da = xr.DataArray(systematic,dims=["times","across_track","channels"],
                          attrs={"long_name":"Systematic uncertainties","units":"Albedo/K"})

    uratio_da = xr.DataArray(irdata['chan_covar_ratio'].values,
                             dims=["times","across_track","ir_channels"],
                             attrs={"long_name":"Channel-to-channel covariance  ratio"})

    #
    # Now merge flags
    #
    uflags = np.zeros(ds['latitude'].shape,dtype=np.int8)
    for i in range(uflags.shape[0]):
        # IR flags are per scanline
        uflags[i,:] = irdata['uncert_flags'].values[i]
        # Vis flags are at pixel level - set 3rd bit for solar contamination
        gd = (visdata['solar_fov_contam'].values[i,:] == 1)
        uflags[i,gd] = (uflags[i,gd]|4)

    uflags_da = xr.DataArray(uflags,dims=["times","across_track"],
                             attrs={"long_name":"Uncertainty flags (bit 1==bad space view (value=1), bit 2==solar contamination of Gain (value=2), bit 3==solar contamination of FOV (value=2)"})

    uncertainties = xr.Dataset(dict(times=time_da,across_track=across_da,
                                    channels=channels_da,
                                    ir_channels=ir_channels_da,
                                    random=random_da,systematic=sys_da,
                                    chan_covar_ratio=uratio_da,
                                    uncert_flags=uflags_da))

    if plot:
        import matplotlib.pyplot as plt
        plt.figure(1)
        plt.subplot(231)
        im=plt.imshow(uncertainties['random'].values[:,:,0])
        plt.title('Chan 1 / Rand')
        plt.colorbar(im)
        plt.subplot(232)
        im=plt.imshow(uncertainties['random'].values[:,:,1])
        plt.title('Chan 2 / Rand')
        plt.colorbar(im)
        plt.subplot(233)
        im=plt.imshow(uncertainties['random'].values[:,:,2])
        plt.title('Chan 3 / Rand')
        plt.colorbar(im)
        plt.subplot(234)
        im=plt.imshow(uncertainties['random'].values[:,:,3])
        plt.title('Chan 4 / Rand')
        plt.colorbar(im)
        plt.subplot(235)
        im=plt.imshow(uncertainties['random'].values[:,:,4])
        plt.title('Chan 5 / Rand')
        plt.colorbar(im)
        if uncertainties['random'].values.shape[2] == 6:
            plt.subplot(236)
            im=plt.imshow(uncertainties['random'].values[:,:,5])
            plt.title('Chan 6 / Rand')
            plt.colorbar(im)
        plt.tight_layout()
            
        plt.figure(2)
        plt.subplot(231)
        im=plt.imshow(uncertainties['systematic'].values[:,:,0])
        plt.title('Chan 1 / Sys')
        plt.colorbar(im)
        plt.subplot(232)
        im=plt.imshow(uncertainties['systematic'].values[:,:,1])
        plt.title('Chan 2 / Sys')
        plt.colorbar(im)
        plt.subplot(233)
        im=plt.imshow(uncertainties['systematic'].values[:,:,2])
        plt.title('Chan 3 / Sys')
        plt.colorbar(im)
        plt.subplot(234)
        im=plt.imshow(uncertainties['systematic'].values[:,:,3])
        plt.title('Chan 4 / Sys')
        plt.colorbar(im)
        plt.subplot(235)
        im=plt.imshow(uncertainties['systematic'].values[:,:,4])
        plt.title('Chan 5 / Sys')
        plt.colorbar(im)
        if uncertainties['random'].values.shape[2] == 6:
            plt.subplot(236)
            im=plt.imshow(uncertainties['systematic'].values[:,:,5])
            plt.title('Chan 6 / Rand')
            plt.colorbar(im)
        plt.tight_layout()
        plt.show()
        
    return uncertainties
        
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('--plot',action='store_true')

    args = parser.parse_args()

    #
    # Read data
    #
    reader_cls = get_reader_class(args.filename)

    reader = reader_cls(tle_dir="/gws/nopw/j04/nceo_uor/users/jmittaz/NPL/AVHRR/TLE",
                        tle_name="TLE_%(satname)s.txt",
                        calibration_method="noaa",
                        adjust_clock_drift=False)
    reader.read(args.filename)
    ds = reader.get_calibrated_dataset()
    mask = reader.mask

    uncert = uncertainty(ds,mask,plot=args.plot)
    print(uncert)
    
