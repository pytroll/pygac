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
import pandas as pd
from importlib.resources import files
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
from pyorbital import astronomy
from pygac import get_reader_class
from pygac.calibration.noaa import Calibrator, calibrate_solar
from pygac.utils import allan_deviation, get_bad_space_counts
from pygac.klm_reader import KLMReader
from pygac.reader import Reader


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


def get_random(noise,av_noise,gain,cal,year,jday,C,D):
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
    Rcal = gain*(C-D)
    #
    # Gain part for all noise sources
    #
    dRcal_dC = gain
    dRcal_dD = -gain

    uncert = (dRcal_dD**2)*(av_noise**2) + (dRcal_dC**2)*(noise**2)

    return np.sqrt(uncert), Rcal

def get_reflectance(Rcal, d_se, sza):
    refl = (Rcal*d_se**2)/np.cos(sza)/100

    return refl

def get_sys(channel, C, D, gain):
    """Get the systematic parts of the vis calibration uncertainty."""
    dRcal_dS = (C-D)
    usys = np.sqrt(0.025**2 + 0.025**2 + 0.01**2 + 0.015**2 + 0.02**2 + 0.025**2)
    # If channel = 2, add water vapour uncertainty
    #
    U_WV = 0.015
    if channel == 2:
        usys_tot = np.sqrt(usys**2 + U_WV**2)
    else:
        usys_tot = usys

    usys_tot *= gain

    uncert = (dRcal_dS**2)*(usys_tot**2)

    return np.sqrt(uncert)

def get_vars(ds,channel):
    """Get variables from xarray"""

    space = ds['vis_space_counts'].values[:,channel]
    counts = ds['channels'].values[:,:,channel-3]

    return space, counts


def vis_uncertainty(ds,mask,plot=False):
    """Create the uncertainty components for the vis channels. These include
    
    1) Random 
        a) Noise
        b) Digitisation
    2) Systematic
        a) MODIS Reflectance uncertainty
        b) SBAF Correction uncertainty
        c) Water Vapour effects
        d) Surface Reflectance variation
        e) Temporal Stability
    
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

    D_1,C_1= get_vars(ds,0)
    if plot:
        plt.subplot(131)
        plt.plot(np.arange(len(D_1)),D_1,',')
        plt.title('0.63$\mu$m')
        plt.ylabel('Space Cnts')
        plt.xlabel('Scanline')

    D_2,C_2 = get_vars(ds,1)
    if plot:
        plt.subplot(132)
        plt.plot(np.arange(len(D_2)),D_2,',')
        plt.title('0.86$\mu$m')
        plt.ylabel('Space Cnts')
        plt.xlabel('Scanline')
        
    if chan_3a:
        D_3,C_3 = get_vars(ds,2)
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
        # Check for bad scanlines and add flag
        #
        if bad_scan[i] == 1:
            print(i, "Bad Space Count Data")
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
        # Get calibration slope
        #
        gain_1 = calibrate_solar(C_1[i, :], 0, year, jday, cal, corr=1)/(C_1[i, :]-D_1[i])
        gain_2 = calibrate_solar(C_2[i, :], 1, year, jday, cal, corr=1)/(C_2[i, :]-D_2[i])
        if chan_3a:
            gain_3 = calibrate_solar(C_3[i, :], 2, year, jday, cal, corr=1)/(C_3[i, :]-D_3[i])

        #
        # Get noise in scaled radiance space
        #
        rad_noise_63, Rcal_1 = get_random(noise1,av_noise1,gain_1,cal,year,jday, C_1[i,:], D_1[i])
        rad_noise_86, Rcal_2 = get_random(noise2,av_noise2,gain_2,cal,year,jday, C_2[i,:], D_2[i])
        if chan_3a:
            rad_noise_12, Rcal_3 = get_random(noise3,av_noise3,gain_3,cal,year,jday, C_3[i,:], D_3[i])


        rcal_rand_63[i,:] = rad_noise_63
        rcal_rand_86[i,:] = rad_noise_86
        if chan_3a:
            rcal_rand_12[i,:] = rad_noise_12

        #
        # Get systematic uncertainty through the measurement equation
        #
        rcal_sys_63[i,:] = get_sys(1, C_1[i,:], D_1[i], gain_1)
        rcal_sys_86[i,:] = get_sys(2, C_2[i,:], D_2[i], gain_2)
        if chan_3a:
            rcal_sys_12[i,:] = get_sys(3, C_3[i,:], D_3[i], gain_3)

    # define flag for solar contamination data
    d_se = ds.attrs["sun_earth_distance_correction_factor"]
    solar_contam_threshold = 0.05
    sza_threshold = 102

    #read = Reader()
    #angles = read.get_angles()
    #sza = angles["sun_zenith"]
    lons, lats = KLMReader()._get_lonlat_from_file()
    times = ds["times"].values
    sza = astronomy.sun_zenith_angle(times[:, np.newaxis],
                                                lons, lats)

    refl_1 = get_reflectance(Rcal_1, d_se, sza)
    refl_2 = get_reflectance(Rcal_2, d_se, sza)
    refl_3 = get_reflectance(Rcal_3, d_se, sza)


    for i in range(len(D_1)):
        contam_pixels_1 = np.where(np.logical_and(refl_1>solar_contam_threshold, sza > sza_threshold),
                                 "Pixel is contaminated", refl_1)

    for i in range(len(D_2)):
        contam_pixels_2 = np.where(np.logical_and(refl_2 > solar_contam_threshold, sza > sza_threshold),
                                 "Pixel is contaminated", refl_2)

    if chan_3a:
        for i in range(len(D_3)):
            contam_pixels_3 = np.where(np.logical_and(refl_3 > solar_contam_threshold, sza > sza_threshold),
                                     "Pixel is contaminated", refl_3)

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

            plt.figure(3)
            plt.subplot(231)
            im = plt.imshow(rcal_rand_63)
            plt.colorbar(im)
            plt.title('0.63$\mu$m (Random)')

            plt.subplot(232)
            im = plt.imshow(rcal_rand_86)
            plt.colorbar(im)
            plt.title('0.86$\mu$m (Random)')

            plt.subplot(233)
            im = plt.imshow(rcal_sys_63)
            plt.colorbar(im)
            plt.title('0.63$\mu$m (Systematic)')

            plt.subplot(234)
            im = plt.imshow(rcal_sys_86)
            plt.colorbar(im)
            plt.title('0.86$\mu$m (Systematic)')

            plt.tight_layout()

        plt.show()

        if chan_3a:
            plt.figure(2)
            plt.subplot(231)
            plt.hist((rcal_rand_63+rcal_sys_63).flatten(), bins=100)
            plt.title('0.63$\mu$m')

            plt.subplot(232)
            plt.hist((rcal_rand_86+rcal_sys_86).flatten(), bins=100)
            plt.title('0.86$\mu$m (Total Uncertainty)')

            plt.subplot(233)
            plt.hist((rcal_rand_12+rcal_sys_12).flatten(), bins=100)
            plt.title('1.2$\mu$m')
            plt.tight_layout()

        else:
            plt.figure(2)
            plt.subplot(221)
            plt.hist((rcal_rand_63+rcal_sys_63).flatten(), bins=100)
            plt.title('0.63$\mu$m (Total Uncertainty)')

            plt.subplot(222)
            plt.hist((rcal_rand_86+rcal_sys_86).flatten(), bins=100)
            plt.title('0.86$\mu$m (Total Uncertainty)')
            plt.tight_layout()

        plt.show()

        if chan_3a:
            cov_chan_rand = np.array([rcal_rand_63[0, :], rcal_rand_86[0, :]], rcal_rand_12[0, :])
            cov = np.cov(cov_chan_rand, bias=True)
            labels = ['0.63$\mu$m', '0.86$\mu$m', '1.2$\mu$m']
            sns.heatmap(cov, annot=True, fmt='g', xticklabels=labels, yticklabels=labels)


        else:
            cov_chan_rand = np.array([rcal_rand_63[0, :], rcal_rand_86[0, :]])
            cov = np.cov(cov_chan_rand, bias=True)
            labels = ['0.63$\mu$m', '0.86$\mu$m']
            sns.heatmap(cov, annot=True, fmt='g', xticklabels=labels, yticklabels=labels)
        plt.title('Covariance Matrix (Random)')
        plt.show()

        if chan_3a:
            cov_chan_sys = np.array([rcal_sys_63[0, :], rcal_sys_86[0, :]], rcal_sys_12[0, :])
            cov = np.cov(cov_chan_sys, bias=True)
            labels = ['0.63$\mu$m', '0.86$\mu$m', '1.2$\mu$m']
            sns.heatmap(cov, annot=True, fmt='g', xticklabels=labels, yticklabels=labels)

        else:
            cov_chan_sys = np.array([rcal_sys_63[0, :], rcal_sys_86[0, :]])
            cov = np.cov(cov_chan_sys, bias=True)
            labels = ['0.63$\mu$m', '0.86$\mu$m']
            sns.heatmap(cov, annot=True, fmt='g', xticklabels=labels, yticklabels=labels, cmap='YlGnBu')
        plt.title('Covariance Matrix (Systematic')
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

    if chan_3a:
        return uncertainties, gain_1, gain_2, gain_3, contam_pixels_1, contam_pixels_2, contam_pixels_3
    else:
        return uncertainties, gain_1, gain_2, contam_pixels_1, contam_pixels_2

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
