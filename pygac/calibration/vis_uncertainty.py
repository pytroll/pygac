#!/usr/bin/env python

# Copyright (c) 2014-2015, 2019 Pytroll Developers

# Author(s):

#   Nicole Yaghnam <nicole.yaghnam@npl.co.uk>

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

"""Uncertainty information based on the noaa.py calibration for VIS channels
"""
from __future__ import division

import numpy as np
import xarray as xr
from importlib.resources import files
import argparse
import matplotlib.pyplot as plt

from pygac import get_reader_class
from pygac.calibration.noaa import Calibrator
from pygac.utils import allan_deviation
from pygac.klm_reader import KLMReader


def get_bad_space_counts(sp_data):
    """Find bad space count data (space count data is voltage clamped so
    should have very close to the same value close to 950 - 960
    Written by J.Mittaz / University of Reading 6 Oct 2024"""

    #
    # Use robust estimators to get thresholds for space counts
    # Use 4 sigma threshold from median
    # Ensure only for good data
    #
    gd = np.isfinite(sp_data)
    if np.sum(gd) == 0:
        sp_bad_data = np.zeros(sp_data.shape,dtype=bool)
        sp_bad_data[:,:] = True
        return sp_bad_data

    quantile = np.quantile(sp_data[gd].flatten(),[0.25,0.75])
    if quantile[0] == quantile[1]:
        quantile[1] = quantile[1]+0.5
    std = (quantile[1]-quantile[0])/1.349
    sp_bad_data = np.zeros(sp_data.shape,dtype=bool)
    sp_bad_data[:,:] = True
    gd = np.isfinite(sp_data)
    sp_bad_data[gd] = ~((np.abs(sp_data[gd] - np.median(sp_data[gd].flatten()))/\
                         std < 5.))

    return sp_bad_data

def get_noise(total_space,window,chan_3a):
    """Get noise estimates from the counts"""
    #
    # Find bad space view data
    #
    bad_data_1 = get_bad_space_counts(total_space[:,:,0])
    bad_data_2 = get_bad_space_counts(total_space[:,:,1])
    if chan_3a:
        bad_data_3 = get_bad_space_counts(total_space[:,:,2])
    bad_scans = np.zeros(total_space.shape[0],dtype=np.int8)
    if chan_3a:
        for i in range(len(bad_scans)):
            if np.any(bad_data_1[i,:]) or np.any(bad_data_2[i,:]) or \
               np.any(bad_data_3[i,:]):
                bad_scans[i] = 1
    else:
        for i in range(len(bad_scans)):
            if np.any(bad_data_1[i,:]) or np.any(bad_data_2[i,:]):
                bad_scans[i] = 1

    #
    # Estimate noise using the Allan deviation plus the digitisation 
    # uncertainty
    #
    # 0.63 micron space counts
    #
    noise1 = allan_deviation(ds["total_vis_space_counts"].values[:, :, 0], bad_scan=bad_scans)
    noise1 = np.sqrt(noise1 * noise1 + 1. / 3)
    #
    # 0.86 micron space counts
    #
    noise2 = allan_deviation(ds["total_vis_space_counts"].values[:, :, 1], bad_scan=bad_scans)
    noise2 = np.sqrt(noise2 * noise2 + 1. / 3)

    #
    # 1.2 micron space counts
    #
    if chan_3a:
        noise3 = allan_deviation(ds["total_vis_space_counts"].values[:, :, 2], bad_scan=bad_scans)
        noise3 = np.sqrt(noise3 * noise3 + 1. / 3)
    else:
        noise3 = None

    #
    # Calculate the uncertainty after averaging - note 10 measurements per
    # scanline
    #
    sqrt_window = np.sqrt(window*10)
    av_noise1 = noise1/sqrt_window
    av_noise2 = noise2/sqrt_window
    if chan_3a:
        av_noise3 = noise3/sqrt_window
    else:
        av_noise3 = None


    return noise1,noise2,noise3,av_noise1,av_noise2,av_noise3,bad_scans


def get_random(noise,av_noise,s0,s1,s2,cal,year,jday,C,D):
    """Get the random parts of the vis calibration uncertainty. Done per
    scanline"""
    #
    # Get time since launch in years
    #
    l_date = Calibrator.date2float(cal.date_of_launch)
    t = (year + jday/365.0) - l_date
    #
    # Measurement Function
    #
    Rcal = s0*(100 + s1*t + s2*t**2)*(C-D)/100
    #
    # Gain part for all noise sources
    #
    dRcal_dC = s0*(100 + s1*t + s2*t**2)/100
    dRcal_dD = -s0*(100 + s1*t + s2*t**2)/100

    uncert = (dRcal_dD**2)*(av_noise**2) + (dRcal_dC**2)*(noise**2)

    return np.sqrt(uncert), Rcal

# def get_sys(channel,
#
#     return np.sqrt(uncert)

def get_vars(ds,channel,wlength,space_threshold,mask):
    """Get variables from xarray"""

    space = ds['vis_space_counts'].values[:,channel]
    counts = ds['channels'].values[:,:,channel-3]

    # Thresholds to flag missing/wrong data for interpolation
    # Remove masked data
    space[mask] = 0
    zeros = space < space_threshold
    nonzeros = np.logical_not(zeros)

    # space[zeros] = np.interp((zeros).nonzero()[0],
    #                          (nonzeros).nonzero()[0],
    #                          space[nonzeros])
    #
    # Make averages and do using pygacs method at this point
    #
    weighting_function = np.ones(wlength, dtype=float) / wlength
    space_convolved = np.convolve(space, weighting_function, "same")

    # take care of the beginning and end
    space_convolved[0:(wlength - 1) // 2] = space_convolved[(wlength - 1) // 2]
    space_convolved[-(wlength - 1) // 2:] = space_convolved[-((wlength + 1) // 2)]

    return space, counts


def vis_uncertainty(ds,mask,plot=False):
    """Create the uncertainty components for the vis channels. These include
    
    1) Random 
          a) Noise
          b) Digitisation
    2) Systematic

    
    Inputs:
          ds : Input xarray dataset containing data for calibration
        mask : pygac mask from reader
    Outputs:
      uncert : xarray dataset containing random and systematic uncertainty
               components
    """
    #
    # Define averaging kernel based on value in noaa.py
    #
    window=51
    space_threshold = 100

    if ds['channels'].values.shape[1] == 409:
        gacdata = True
    else:
        gacdata = False

    avhrr_name = ds.attrs["spacecraft_name"]

    #
    # Test for channel 3a
    #

    if KLMReader._get_vis_channels_to_calibrate == [0,1,2]:
        if ds['channels'].value == np.nan:
            chan_3a = False
        else:
            chan_3a = True
    else:
        chan_3a = False

    #
    # Get calibration coefficients
    #
    cal = Calibrator(
        ds.attrs["spacecraft_name"])
    s0_1 = cal.s0[0]
    s1_1 = cal.s0[1]
    s2_1 = cal.s0[2]
    s0_2 = cal.s1[0]
    s1_2 = cal.s1[1]
    s2_2 = cal.s1[2]
    s0_3 = cal.s2[0]
    s1_3 = cal.s2[1]
    s2_3 = cal.s2[2]

    #
    # Get time since launch in years
    #
    times = ds.coords["times"]
    start_time = times[0].dt
    year = start_time.year.item()
    jday = start_time.dayofyear.item()

    #
    # Get variables for 10 sampled case
    #
    total_space = ds['total_vis_space_counts'].values[:,:,:]
    
    #
    # Noise elements
    #
    noise1,noise2,noise3,av_noise1,av_noise2,av_noise3,bad_scan \
               = get_noise(total_space,window,chan_3a)

    #
    # Get variables used on the calibration
    #
    if plot:
        plt.figure(1)

    D_1,C_1= get_vars(ds,0,window,space_threshold,mask)
    if plot:
        plt.subplot(131)
        plt.plot(np.arange(len(D_1)),D_1,',')
        plt.title('0.63$\mu$m')
        plt.ylabel('Space Cnts')
        plt.xlabel('Scanline')

    D_2,C_2 = get_vars(ds,1,window,space_threshold,mask)
    if plot:
        plt.subplot(132)
        plt.plot(np.arange(len(D_2)),D_2,',')
        plt.title('0.86$\mu$m')
        plt.ylabel('Space Cnts')
        plt.xlabel('Scanline')
        
    if chan_3a:
        D_3,C_3 = get_vars(ds,2,window,space_threshold,mask)
        if plot:
            plt.subplot(133)
            plt.plot(np.arange(len(D_3)),D_3,',')
            plt.title('1.2$\mu$m')
            plt.ylabel('Space Cnts')
            plt.xlabel('Scanline')
    if plot:
        plt.tight_layout()

    #
    # Systematic components
    #

    #
    # Loop round scanlines
    #
    rcal_rand_63 = np.zeros(C_2.shape,dtype=C_2.dtype)
    rcal_rand_86 = np.zeros(C_2.shape,dtype=C_2.dtype)
    rcal_rand_12 = np.zeros(C_2.shape,dtype=C_2.dtype)
    if not chan_3a:
        rcal_rand_12[:,:] = np.nan
    rcal_sys_63 = np.zeros(C_2.shape,dtype=C_2.dtype)
    rcal_sys_86 = np.zeros(C_2.shape,dtype=C_2.dtype)
    rcal_sys_12 = np.zeros(C_2.shape,dtype=C_2.dtype)
    if not chan_3a:
        rcal_sys_12[:,:] = np.nan
    for i in range(len(D_2)):
        #
        # Check for bad scanlines
        #
        if bad_scan[i] == 1:
            rcal_rand_63[i,:] = np.nan
            rcal_rand_86[i,:] = np.nan
            if chan_3a:
                rcal_rand_12[i,:] = np.nan
            rcal_sys_63[i,:] = np.nan
            rcal_sys_86[i,:] = np.nan
            if chan_3a:
                rcal_sys_12[i,:] = np.nan
            continue



        #
        # Get noise in scaled radiance space  (rename)
        #
        rad_noise_63, Rcal_1 = get_random(noise1,av_noise1,s0_1,s1_1,s2_1,cal,year,jday, C_1[i,:], D_1[i])
        rad_noise_86, Rcal_2 = get_random(noise2,av_noise2,s0_2,s1_2,s2_2,cal,year,jday, C_2[i,:], D_2[i])
        if chan_3a:
            rad_noise_12, Rcal_3 = get_random(noise3,av_noise3,s0_3,s1_3,s2_3,cal,year,jday, C_3[i,:], D_3[i])


        rcal_rand_63[i,:] = rad_noise_63
        rcal_rand_86[i,:] = rad_noise_86
        if chan_3a:
            rcal_rand_12[i,:] = rad_noise_12

        #
        # Get systematic uncertainty through the measurement equation
        #

    if plot:
        if chan_3a:
            plt.figure(2)
            plt.subplot(231)
            plt.hist(rcal_rand_63.flatten(),bins=100)
            plt.title('0.63$\mu$m')
            
            plt.subplot(232)
            plt.hist(rcal_rand_86.flatten(),bins=100)
            plt.title('0.86$\mu$m (Random)')

            plt.subplot(233)
            plt.hist(rcal_rand_12.flatten(),bins=100)
            plt.title('1.2$\mu$m')
            
            plt.subplot(234)
            plt.hist(rcal_sys_63.flatten(),bins=100)
            plt.title('0.63$\mu$m')
            
            plt.subplot(235)
            plt.hist(rcal_sys_86.flatten(),bins=100)
            plt.title('0.86$\mu$m (Systematic)')
            plt.xlabel('Uncertainty')

            plt.subplot(236)
            plt.hist(rcal_sys_12.flatten(),bins=100)
            plt.title('1.2$\mu$m')
            plt.tight_layout()

            plt.figure(3)
            plt.subplot(231)
            im=plt.imshow(rcal_rand_63)
            plt.colorbar(im)
            plt.title('0.63$\mu$m')

            plt.subplot(232)
            im=plt.imshow(rcal_rand_86)
            plt.colorbar(im)
            plt.title('0.86$\mu$m (Random)')

            plt.subplot(233)
            im=plt.imshow(rcal_rand_12)
            plt.colorbar(im)
            plt.title('1.2$\mu$m')

            plt.subplot(234)
            im=plt.imshow(rcal_sys_63)
            plt.colorbar(im)
            plt.title('0.63$\mu$m')

            plt.subplot(235)
            im=plt.imshow(rcal_sys_86)
            plt.colorbar(im)
            plt.title('0.86$\mu$m (Systematic)')

            plt.subplot(236)
            im=plt.imshow(rcal_sys_12)
            plt.colorbar(im)
            plt.title('1.2$\mu$m')
            plt.tight_layout()
        else:
            plt.figure(2)
            plt.subplot(221)
            plt.hist(rcal_rand_63.flatten(),bins=100)
            plt.title('0.63$\mu$m (Random)')
            
            plt.subplot(222)
            plt.hist(rcal_rand_86.flatten(),bins=100)
            plt.title('0.86$\mu$m (Random)')

            plt.subplot(223)
            plt.hist(rcal_sys_63.flatten(),bins=100)
            plt.title('0.63$\mu$m (Systematic)')
            plt.xlabel('Uncertainty')
            
            plt.subplot(224)
            plt.hist(rcal_sys_86.flatten(),bins=100)
            plt.title('0.86$\mu$m (Systematic)')
            plt.xlabel('Uncertainty')

            plt.tight_layout()
        plt.show()
    #
    # Output uncertainties
    #
    random = np.zeros((rcal_rand_86.shape[0],rcal_rand_86.shape[1],3))
    systematic = np.zeros((rcal_rand_86.shape[0],rcal_rand_86.shape[1],3))

    random[:,:,0] = rcal_rand_63
    random[:,:,1] = rcal_rand_86
    if chan_3a:
        random[:,:,2] = rcal_rand_12
    else:
        random[:,:,2] = np.nan
    systematic[:,:,0] = rcal_sys_63
    systematic[:,:,1] = rcal_sys_86
    if chan_3a:
        systematic[:,:,2] = rcal_sys_12
    else:
        systematic[:,:,2] = np.nan


    time = (ds["times"].values - np.datetime64("1970-01-01 00:00:00"))/\
           np.timedelta64(1,'s')
    time_da = xr.DataArray(time,dims=["times"],attrs={"long_name":"scanline time",\
                                                     "units":"seconds since 1970-01-01"})
    across_da = xr.DataArray(np.arange(random.shape[1]),dims=["across_track"])
    vis_channels_da = xr.DataArray(np.array([1,2,3]),dims=["vis_channels"])
    random_da = xr.DataArray(random,dims=["times","across_track","vis_channels"],\
                             attrs={"long_name":"Random uncertainties","units":""})
    sys_da = xr.DataArray(systematic,dims=["times","across_track","vis_channels"],\
                          attrs={"long_name":"Systematic uncertainties","units":""})

    uncertainties = xr.Dataset(dict(times=time_da,across_track=across_da,\
                               vis_channels=vis_channels_da,\
                                    random=random_da,systematic=sys_da))

    return uncertainties


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('filename')

    args = parser.parse_args()

    #
    # Read data
    #
    reader_cls = get_reader_class(args.filename)

    reader = reader_cls(tle_dir="/gws/nopw/j04/npl_eo/users/nyaghnam/pygac/gapfilled_tles",
                        tle_name="TLE_%(satname).txt",
                        calibration_method="noaa",
                        adjust_clock_drift=False)
    reader.read(args.filename)
    ds = reader.get_calibrated_dataset()
    mask = reader.mask
    uncert = vis_uncertainty(ds,mask,plot=True)
    print(uncert)
