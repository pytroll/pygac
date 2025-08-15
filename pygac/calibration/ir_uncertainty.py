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
from importlib.resources import files
import argparse
from scipy.optimize import curve_fit

from pygac import get_reader_class
from pygac.calibration.noaa import Calibrator
from pygac.calibration.noaa import get_prt_nos

def allan_deviation(space,bad_scan=None):
    """Determine the Allan deviation (noise) from space view counts filtering
    out bad space view lines. Written by J.Mittaz, University of Reading"""

    if len(space.shape) != 2:
        raise Exception("utils.allan_deviation input space view not 2-dimensional")
    #
    # Get good scanlines if filter present
    #
    if bad_scan is not None:
        gd = (bad_scan == 0)
        newsp = space[gd,:]
    else:
        newsp = space

    #
    # Allan deviation is sqrt of allan variance which is
    #
    #        allan_variance = 0.5 <(y_n+1-y_n)**2>
    #
    allan = 0.
    nline = newsp.shape[0]
    nscan = newsp.shape[1]
    #
    # Sum of squared pixel differences
    #
    for i in range(nline):
        allan += np.sum((newsp[i,1:]-newsp[i,0:-1])**2)
        
    #
    # Mean of squared pixel to pixel difference (expectation)
    #
    allan /= (nline*(nscan-1))
    #
    # And 1/2 term
    #
    allan *= 0.5

    #
    # Return the Allan deviation in counts (sqrt Allan variance)
    #
    return np.sqrt(allan)

class convBT(object):
    """Routine to covert temperature to radiance and visa-versa"""
    def t_to_rad(self,tprt):

        tsBB = self.A + self.B*tprt
        return self.nBB_num / (np.exp(self.c2_nu_c / tsBB) - 1.0)
        
    def rad_to_t(self,rad):

        corrT = self.c2_nu_c/np.log((self.nBB_num/rad)+1.)
        return (corrT-self.A)/self.B

    def rad_to_t_uncert(self,rad,urad):

        T = self.rad_to_t(rad)
        T1 = self.rad_to_t(rad+urad)
        T2 = self.rad_to_t(rad-urad)

        return T,(np.abs(T-T1)+np.abs(T-T2))/2.

    def t_to_rad_uncert(self,T,uT):

        rad = self.t_to_rad(T)
        rad1 = self.t_to_rad(T+uT)
        rad2 = self.t_to_rad(T-uT)

        return rad,(np.abs(rad-rad1)+np.abs(rad-rad2))/2.

    def __init__(self,cal,chan):

        # constants
        self.c1 = 1.1910427e-5  # mW/m^2/sr/cm^{-4}
        self.c2 = 1.4387752  # cm K
        # coefficients
        self.A = cal.to_eff_blackbody_intercept[chan]
        self.B = cal.to_eff_blackbody_slope[chan]
        self.nu_c = cal.centroid_wavenumber[chan]
        self.nBB_num = self.c1 * (self.nu_c**3)
        self.c2_nu_c = self.c2 * self.nu_c

def get_bad_space_counts(sp_data,ict_data=None):
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
    if ict_data is not None:
        sp_bad_data[gd] = ~((np.abs(sp_data[gd] - np.median(sp_data[gd].flatten()))/\
                             std < 5.)&(ict_data[gd] > 0))
    else:
        sp_bad_data[gd] = ~((np.abs(sp_data[gd] - np.median(sp_data[gd].flatten()))/\
                             std < 5.))

    return sp_bad_data

def get_noise(total_space,total_ict,window,twelve_micron):
    """Get noise estimates from the counts"""

    #
    # Find bad space view data
    #
    bad_data_1 = get_bad_space_counts(total_space[:,:,0],
                                      ict_data=total_ict[:,:,0])
    bad_data_2 = get_bad_space_counts(total_space[:,:,1],
                                      ict_data=total_ict[:,:,1])
    if twelve_micron:
        bad_data_3 = get_bad_space_counts(total_space[:,:,2],
                                          ict_data=total_ict[:,:,2])
    bad_scans = np.zeros(total_space.shape[0],dtype=np.int8)
    if twelve_micron:
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
    # 3.7 micron space counts
    #
    noise1 = allan_deviation(total_space[:,:,0],bad_scan=bad_scans)
    noise1 = np.sqrt(noise1*noise1 + 1./3)

    #
    # 11 micron space counts
    #
    noise2 = allan_deviation(total_space[:,:,1],bad_scan=bad_scans)
    noise2 = np.sqrt(noise2*noise2 + 1./3)

    #
    # 12 micron space counts
    #
    if twelve_micron:
        noise3 = allan_deviation(total_space[:,:,2],bad_scan=bad_scans)
        noise3 = np.sqrt(noise3*noise3 + 1./3)
    else:
        noise3 = None

    #
    # 3.7 micron ICT counts
    #
    ict_noise1 = allan_deviation(total_ict[:,:,0],bad_scan=bad_scans)
    ict_noise1 = np.sqrt(ict_noise1*ict_noise1 + 1./3)

    #
    # 11 micron ICT counts
    #
    ict_noise2 = allan_deviation(total_ict[:,:,1],bad_scan=bad_scans)
    ict_noise2 = np.sqrt(ict_noise2*ict_noise2 + 1./3)

    #
    # 12 micron ICT counts
    #
    if twelve_micron:
        ict_noise3 = allan_deviation(total_ict[:,:,2],bad_scan=bad_scans)
        ict_noise3 = np.sqrt(ict_noise3*ict_noise3 + 1./3)
    else:
        ict_noise3 = None

    #
    # Calculate the uncertainty after averaging - note 10 measurements per
    # scanline
    #
    sqrt_window = np.sqrt(window*10)
    av_noise1 = noise1/sqrt_window
    av_noise2 = noise2/sqrt_window
    if twelve_micron:
        av_noise3 = noise3/sqrt_window
    else:
        av_noise3 = None
    av_ict_noise1 = ict_noise1/sqrt_window
    av_ict_noise2 = ict_noise2/sqrt_window
    if twelve_micron:
        av_ict_noise3 = ict_noise3/sqrt_window
    else:
        av_ict_noise3 = None

    return noise1,noise2,noise3,av_noise1,av_noise2,av_noise3,\
        av_ict_noise1,av_ict_noise2,av_ict_noise3,bad_scans

def smooth_data(y,length):
    """Smooth data over given length"""
    
    outy = np.zeros(len(y),dtype=y.dtype)
    leny = len(y)-1
    for i in range(leny+1):
        minx=max([i-length,0])
        maxx=min([i+length,leny])
        outy[i] = np.mean(y[minx:maxx+1])

    return outy

def get_uICT(gainval,CS,CICT,Tict,NS,convT,bad_scans,solar_scans,window,
             plt=None,nfigure=None):
    """Get ICT temperature gradient uncertainty based on analysis of the
    gain variations in the 3.7 micron channel
    Estimate on "window" scanline length
    Only uses 'good' data"""

    #
    # Only use good data
    #
    gd = (bad_scans == 0)&(solar_scans == 0)
    Lict = convT.t_to_rad(Tict[gd])

    gain = (Lict-NS)/(CS[gd]-CICT[gd])
    sp_ict = CS[gd]-CICT[gd]
    Tcorr = convT.A+convT.B*Tict[gd]
    
    dGain_dICT = (convT.B/sp_ict)*\
        convT.nBB_num*np.exp(convT.c2_nu_c/Tcorr)*\
        (convT.c2_nu_c/(Tcorr**2))/\
        (np.exp(convT.c2_nu_c/Tcorr)-1.)**2
    
    dgain = (gain-gainval)

    dT = np.zeros(len(CS))
    dT[:] = np.nan
    #
    # Calculate delta ICT for 'good' cases
    #
    dT[gd] = dgain/dGain_dICT

    uICT = np.zeros(len(CS))    
    for i in range(len(dT)):
        if np.isfinite(dT[i]):
            minx = max([0,i-window//2])
            maxx = min([len(dT)+1,i+window//2+1])
            deltaT = dT[minx:maxx]
            gd = np.isfinite(deltaT)
            if np.sum(gd) > 0:
                uICT[i] = np.std(deltaT[gd])
            else:
                uICT[i] = np.nan
        else:
            uICT[i] = np.nan

    if plt is not None:
        nfigure = nfigure+1
        plt.figure(nfigure)
        plt.subplot(311)
        plt.plot(np.arange(len(gain)),gain*1e3,'.')
        plt.axhline(y=gainval*1e3,linestyle='--',color='red')
        plt.ylabel('3.7$\mu$m Gain $\\times 10^{-3}$')
        plt.subplot(312)
        plt.plot(np.arange(len(dT)),dT,'.')
        plt.ylabel('$\Delta$T / K')
        plt.subplot(313)
        plt.plot(np.arange(len(uICT)),uICT,'.')
        plt.xlabel('Scanline')
        plt.ylabel('uICT / K')
        plt.tight_layout()
        
    return uICT,nfigure

def get_ict_uncert(Tict,prt_random,prt_bias,uICT,convT):
    """Get uncertainty in radiance/temperature of ICT on the basis of 
    ICT uncertainties"""

    #
    # Note in operational calibration prts are averaged over 4
    # so input prt uncertainties need to be divided by 2
    #
    rad,urand = convT.t_to_rad_uncert(Tict,prt_random/2.)

    #
    # Systematic doesn't average down
    #
    prt_sys = np.sqrt(prt_bias**2 + uICT**2)

    return urand,prt_sys,rad

class fit_ict_pars(object):
    '''Model for ICT gradients'''
    
    def model(self,X,a0,a1,a2,a3):
        
        ICT = X[:,0] + a0*X[:,1] + a1*X[:,2] + a2*X[:,3] + a3*X[:,4]
        rad_ict = self.convT.t_to_rad(ICT)

        return (rad_ict - X[:,7])/(X[:,5] - X[:,6])
    
    def __init__(self,convT):

        self.convT = convT

def find_solar_ind(position,new_gain2,new_cs,new_ct,new_ict,new_ict1,
                   new_ict2,new_ict3,new_ict4,solZA,first=True,
                   plot=False):
    '''Find possible solar contamination on one side of the solar max
    position'''
    
    pos = np.nonzero(solZA == solZA.max())[0][0]
    #
    # Look before or after max solZA location
    #
    if first:
        posmin = np.nonzero(solZA[0:pos] == solZA[0:pos].min())[0][0]
    else:
        posmin = pos+np.nonzero(solZA[pos:] == solZA[pos:].min())[0][0]
    if pos < posmin:
        posn = position[pos:posmin]
        newg = new_gain2[pos:posmin]
        new_cs2 = new_cs[pos:posmin] 
        new_ct2 = new_ct[pos:posmin] 
        new_ict_2 = new_ict[pos:posmin] 
        new_ict2_1 = new_ict1[pos:posmin] 
        new_ict2_2 = new_ict2[pos:posmin]
        new_ict2_3 = new_ict3[pos:posmin] 
        new_ict2_4 = new_ict4[pos:posmin] 
        sza = solZA[pos:posmin]
    else:
        posn = position[posmin:pos]
        newg = new_gain2[posmin:pos]
        new_cs2 = new_cs[posmin:pos] 
        new_ct2 = new_ct[posmin:pos] 
        new_ict_2 = new_ict[posmin:pos] 
        new_ict2_1 = new_ict1[posmin:pos] 
        new_ict2_2 = new_ict2[posmin:pos] 
        new_ict2_3 = new_ict3[posmin:pos] 
        new_ict2_4 = new_ict4[posmin:pos] 
        sza = solZA[posmin:pos]

    #
    # Only look in range satellite ZA >= 100 and <= 125
    #
    gd = (sza >= 100)&(sza <= 125)&np.isfinite(newg)
    if np.sum(gd) == 0:
        print('ERROR: No data available to find solar contamination')
        print('       MinSZA = {0:f} MaxSZA = {1:f}'.format(sza.min(),sza.max()))
        return None,None,None
    posn2 = posn[gd]
    pos2 = np.nonzero(newg[gd] == newg[gd].min())[0][0]
    newg = newg[gd]
    peak_gain = newg[pos2]
    #
    # Track back to a defined limit when gradient is reversed or zero
    #
    peak_location = posn2[pos2]
    ok1 = False
    startpos = 0
    endpos = len(newg)-1
    window_size = 40
    while not ok1:
        side1 = 0
        for i in range(pos2,0,-1):
            if i > pos2-window_size:
                continue
            minx = max([0,i-window_size])
            maxx = min([peak_location,i+window_size])
            y = newg[minx:maxx+1]
            x = np.arange(len(y))
            p = np.polyfit(x,y,1)
            if p[0] >= 0.:
                side1 = i
                break
        #
        # Check that we are not at an edge
        #
        if side1 == 0:
            if plot:
                plt.show()
            print('Cannot find peak for solar contamination')
            return None,None,None
        
        if pos2+window_size < len(newg):
            ok1 = True
        else:
            #
            # Find next peak < side1 position
            #
            pos2 = np.nonzero(newg[:side1] == newg[:side1].min())[0][0]
            peak_location = posn2[pos2]

    #
    # Look at other side - shouldn't have same problem as side1 with possible
    # erroneous peak
    #
    side2 = len(newg)-1
    ok2 = False
    for i in range(pos2,len(newg)):
        if i < pos2+window_size:
            continue
        minx = max([0,i-window_size])
        maxx = min([peak_location,i+window_size])
        y = newg[minx:maxx+1]
        x = np.arange(len(y))
        p = np.polyfit(x,y,1)
        if p[0] <= 0.:
            ok2 = True
            side2 = i
            break

    #
    # Add Take larger difference case
    #
    diff1 = np.abs(side1-pos2)
    diff2 = np.abs(side2-pos2)
    if diff1 > diff2:
        side1 = diff1
        side2 = diff1
    else:
        side1 = diff2
        side2 = diff2

    #
    # Check to see if we are already at a background value or have to go
    # further (case where second peak on side of signal
    #
    newpos1 = max([0,pos2-side1])
    newpos2 = min([len(newg)-1,pos2+side2])
    maxdiff1 = np.abs(peak_gain-newg[newpos1])
    maxdiff2 = np.abs(peak_gain-newg[newpos2])
    maxdiff = max([maxdiff1,maxdiff2])
    update_side1 = False
    for i in range(side1+side1//2):
        bloc = pos2-side1-i
        if bloc <= 0:
            break
        test_ratio = maxdiff/np.abs(peak_gain-newg[bloc])
        if test_ratio > 0.7 and test_ratio < 1.3:
            update_side1 = True
            side1 = bloc
            break

    update_side2 = False
    for i in range(side2+side2//2):
        bloc = pos2+side2+i
        if bloc >= len(newg):
            break
        test_ratio = maxdiff/np.abs(peak_gain-newg[bloc])
        if test_ratio > 0.7 and test_ratio < 1.3:
            update_side2 = True
            side2 = bloc
            break

    #
    # Reset to full location
    #
    if not update_side1:
        side1 = newpos1
    if not update_side2:
        side2 = newpos2
    
    #
    # Can't get background so eject
    #
    if side1 == 0:
        print('ERROR: Cannot get a background for solar contamination')
        return None,None,None

    #
    # Background and sigma
    #
    maxx = side1
    minx = max([0,maxx-300])
    backdata1 = newg[minx:maxx]
    pos_1 = (minx+maxx)/2.
    gd = np.isfinite(backdata1)
    if np.sum(gd) == 0:
        print('ERROR: Cannot get background for solar contamination')
        return None,None,None
    backg1 = np.mean(backdata1[gd])
    backg1_std = np.std(backdata1[gd])

    maxx = side2
    minx = min([len(newg)-1,maxx+300])
    backdata2 = newg[maxx:minx]
    gd = np.isfinite(backdata2)
    if np.sum(gd) == 0:
        pos_2 = -1
    else:
        backg2 = np.mean(backdata2[gd])
        pos_2 = (minx+maxx)/2.
        backg2_std = np.std(backdata2[gd])    
        
    #
    # If pos_2 there, and OK interpolate
    #
    if pos_2 > -1:
        #
        # Make sure back2 is < 50% of peak
        #
        if (backg1-peak_gain)*0.25 > (backg2-peak_gain):
            #
            # Less than 25% from peak so dont continue
            #
            print('ERROR: back2 is too high (solar contamination)')
            return None,None,None
        slope = (backg2-backg1)/(pos_2-pos_1)
        cnst = backg1 - slope*pos_1
        backg = cnst + slope*pos2
        backg_std = (backg1_std+backg2_std)/2.
    else:
        print('ERROR: No back2 value (solar contamination)')
        return None,None,None

    sigma = np.abs(peak_gain-backg)/backg_std
    if sigma <= 7.5:
        return None,None,None
    
    side1 = posn2[side1]
    side2 = posn2[side2]

    return side1,side2,peak_location
        
def find_solar(ds,mask,convT=None,outgain=False,plot=False,out_time=False):
    '''Find solar contamination/variable gain points for GAC and for ~full
    orbits'''

    #
    # Get parameters for kernals/uncertainties
    #
    window, prt_bias, prt_sys, prt_threshold, ict_threshold, \
        space_threshold = get_uncert_parameter_thresholds()
    if ds['channels'].values.shape[1] == 409:
        gacdata = True
        scale = 1
    else:
        gacdata = False
        scale = 3
        #
        # Code only works for GAC
        #
        print('ERROR: Solar contamination detection only works for GAC data')
        if outgain:
            return None,None,None,None,None,None
        else:
            return None,None,None,None
        
    avhrr_name = ds.attrs["spacecraft_name"]
    #
    # Get calibration coefficients for 3.7 micron channel
    #
    cal = Calibrator(
        ds.attrs["spacecraft_name"])
    NS_1 = cal.space_radiance[0]
    c0_1 = cal.b[0,0]
    c1_1 = cal.b[0,1]
    c2_1 = cal.b[0,2]

    if convT is None:
        convT = convBT(cal,0)

    if out_time:
        CS_1,CICT_1,CE_1,Tict,ict1,ict2,ict3,ict4,solZAin,time \
            = get_vars(ds,0,convT,
                       window,
                       prt_threshold,
                       ict_threshold,
                       space_threshold,
                       gacdata,
                       cal,
                       mask,
                       out_prt=True,
                       out_solza=True,
                       out_time=True)
    else:
        CS_1,CICT_1,CE_1,Tict,ict1,ict2,ict3,ict4,solZAin \
            = get_vars(ds,0,convT,
                       window,
                       prt_threshold,
                       ict_threshold,
                       space_threshold,
                       gacdata,
                       cal,
                       mask,
                       out_prt=True,
                       out_solza=True,
                       out_time=False)

    meanT = (ict1+ict2+ict3+ict4)/4.
    radBB = convT.t_to_rad(meanT)
    gain3 = (radBB-NS_1)/(CS_1-CICT_1)

    #
    # Get location of min PRT stdev as closest to zero error case
    #
    X = np.zeros((len(meanT),4))
    X[:,0] = ict1
    X[:,1] = ict2
    X[:,2] = ict3
    X[:,3] = ict4
    stdev = np.std(X,axis=1)
    try:
        minstd_pos = np.nonzero(stdev == stdev.min())[0][0]
    except:
        minstd_pos = np.nonzero(stdev == stdev.min())[0]

    #
    # Get PRT differences from Mean
    #
    prt_diff1 = ict1 - meanT
    prt_diff2 = ict2 - meanT
    prt_diff3 = ict3 - meanT
    prt_diff4 = ict4 - meanT

    #
    # Now fit model
    #
    X = np.zeros((len(meanT),8))
    Y = np.zeros((len(meanT)))
    X[:,0] = meanT[:]
    X[:,1] = prt_diff1[:]
    X[:,2] = prt_diff2[:]
    X[:,3] = prt_diff3[:]
    X[:,4] = prt_diff4[:]
    X[:,5] = CS_1[:]
    X[:,6] = CICT_1[:]
    X[:,7] = NS_1
    Y[:] = gain3[minstd_pos]

    model = fit_ict_pars(convT)
    p,covar = curve_fit(model.model,X,Y,p0=[0.,0.,0.,0.])
    
    #
    # Get updated gain
    #
    new_gain = model.model(X,p[0],p[1],p[2],p[3])

    #
    # Find peak gain within solZA 100-125 degrees
    #
    # First find range between min/max values
    #
    position = np.arange(len(solZAin)).astype(dtype=np.int32)
    gd = np.isfinite(solZAin)
    position = position[gd]
    new_gain2 = new_gain[gd]
    new_cs = CS_1[gd]
    new_ct = CICT_1[gd]
    new_ict = Tict[gd]
    new_ict1 = ict1[gd]
    new_ict2 = ict2[gd]
    new_ict3 = ict3[gd]
    new_ict4 = ict4[gd]                    
    gain3_2 = gain3[gd]
    solZA = solZAin[gd]

    #
    # Now look before and after solZA max position
    #
    #                                                                          
    # One side of solZA max
    #                                                      
    side1_1,side2_1,peak_location_1 = find_solar_ind(position,new_gain2,
                                                     new_cs,new_ct,new_ict,
                                                     new_ict1,new_ict2,
                                                     new_ict3,new_ict4,
                                                     solZA,first=True,
                                                     plot=plot)
    # 
    # Other side of solZA max
    #
    side1_2,side2_2,peak_location_2 = find_solar_ind(position,new_gain2,
                                                     new_cs,new_ct,new_ict,
                                                     new_ict1,new_ict2,
                                                     new_ict3,new_ict4,
                                                     solZA,first=False,
                                                     plot=plot)
    
    #
    # Plot PRT variation and gain differences (before/after)
    #
    if plot:
        plt.figure(1)
        x = np.arange(len(prt_diff1))
        plt.subplot(211)
        plt.plot(x,prt_diff1,label='PRT1')
        plt.plot(x,prt_diff2,label='PRT2')
        plt.plot(x,prt_diff3,label='PRT3')
        plt.plot(x,prt_diff4,label='PRT4')
        plt.ylabel('Temp. Diff / K')
        plt.legend()
        plt.subplot(212)
        plt.plot(x,gain3)
        plt.ylabel('3.7$\mu$m Gain $\\times 10^{-3}$')
        plt.axhline(y=gain3[minstd_pos],linestyle='--',color='red')
        plt.figure(2)
        x = np.arange(len(prt_diff1))
        plt.subplot(211)
        plt.plot(x,prt_diff1,label='PRT1')
        plt.plot(x,prt_diff2,label='PRT2')
        plt.plot(x,prt_diff3,label='PRT3')
        plt.plot(x,prt_diff4,label='PRT4')
        plt.xlabel('Scanline')
        plt.ylabel('Temp. Diff / K')
        plt.legend()
        plt.subplot(212)
        plt.plot(x,gain3*1e3,label='Orig')
        plt.plot(x,new_gain*1e3,label='New')
        plt.xlabel('Scanline')
        plt.ylabel('3.7$\mu$m Gain $\\times 10^{-3}$')
        plt.legend()
        if side1_1 is not None:
            plt.axvline(x=side1_1,color='red')
            plt.axvline(x=side2_1,color='red')
        if side1_2 is not None:
            plt.axvline(x=side1_2,color='red')
            plt.axvline(x=side2_2,color='red')
        plt.tight_layout()
        plt.figure(3)
        plt.plot(x,gain3*1e3)
        plt.xlabel('Scanline')
        plt.ylabel('3.7$\mu$m Gain $\\times 10^{-3}$')
        if side1_1 is not None:
            plt.axvline(x=side1_1,color='red')
            plt.axvline(x=side2_1,color='red')
        if side1_2 is not None:
            plt.axvline(x=side1_2,color='red')
            plt.axvline(x=side2_2,color='red')
        plt.show()
        
    if outgain:
        if side1_1 is not None and side1_2 is not None:
            if out_time:
                return time[side1_1],time[side2_1],time[peak_location_1],\
                    solZAin[peak_location_1],time[side1_2],time[side2_2],\
                    time[peak_location_2],solZAin[peak_location_2],\
                    time[minstd_pos],gain3[minstd_pos]
            else:
                return side1_1,side2_1,peak_location_1,\
                    solZAin[peak_location_1],side1_2,side2_2,\
                    peak_location_2,solZAin[peak_location_2],\
                    minstd_pos,gain3[minstd_pos]
        elif side1_1 is not None and side1_2 is None:
            if out_time:
                return time[side1_1],time[side2_1],time[peak_location_1],\
                    solZAin[peak_location_1],None,None,None,None,\
                    time[minstd_pos],gain3[minstd_pos]
            else:
                return side1_1,side2_1,peak_location_1,\
                    solZAin[peak_location_1],None,None,None,None,\
                    minstd_pos,gain3[minstd_pos]
        elif side1_1 is None and side1_2 is not None:
            if out_time:
                return None,None,None,None,time[side1_2],time[side2_2],\
                    time[peak_location_2],solZAin[peak_location_2],\
                    time[minstd_pos],gain3[minstd_pos]
            else:
                return None,None,None,None,side1_2,side2_2,\
                    peak_location_2,solZAin[peak_location_2],\
                    minstd_pos,gain3[minstd_pos]
        else:
            if out_time:
                return None,None,None,None,None,None,None,None,\
                    time[minstd_pos],gain3[minstd_pos]
            else:
                return None,None,None,None,None,None,None,None,\
                    minstd_pos,gain3[minstd_pos]
    else:
        if side1_1 is not None and side1_2 is not None:
            if out_time:
                return time[side1_1],time[side2_1],time[peak_location_1],\
                    solZAin[peak_location_1],time[side1_2],time[side2_2],\
                    time[peak_location_2],solZAin[peak_location_2]
            else:
                return side1_1,side2_1,peak_location_1,\
                    solZAin[peak_location_1],side1_2,side2_2,\
                    peak_location_2,solZAin[peak_location_2]
        elif side1_1 is not None and side1_2 is None:
            if out_time:
                return time[side1_1],time[side2_1],time[peak_location_1],\
                    solZAin[peak_location_1],None,None,None,None
            else:
                return side1_1,side2_1,peak_location_1,\
                    solZAin[peak_location_1],None,None,None,None
        elif side1_1 is None and side1_2 is not None:
            if out_time:
                return None,None,None,None,time[side1_2],time[side2_2],\
                    time[peak_location_2],solZAin[peak_location_2]
            else:
                return None,None,None,None,side1_2,side2_2,\
                    peak_location_2,solZAin[peak_location_2]
        else:            
            return None,None,None,None,None,None,None,None
        
def get_random(channel,noise,av_noise,ict_noise,ict_random,Lict,CS,CE,CICT,NS,
               c1,c2):
    """Get the random parts of the IR calibration uncertainty. Done per 
    scanline"""
    #
    # Gain part for all noise sources
    #
    dLlin_dCS = (Lict-NS)*(CS-CE)/(CS-CICT)**2 + (Lict-NS)/(CS-CICT)
    dLlin_dCICT = -(Lict-NS)*(CS-CE)/(CS-CICT)**2 
    dLlin_dCE = -(Lict-NS)/(CS-CICT)
    dLlin_dLict = (CS-CE)/(CS-CICT)

    #
    # If channel = 2,3 (11/12) then add non-linear part
    #
    if channel == 1:
        uncert = (dLlin_dCS**2)*(av_noise**2) + \
                 (dLlin_dCICT**2)*(ict_noise**2) + \
                 (dLlin_dCE**2)*(noise**2) + \
                 (dLlin_dLict**2)*(ict_random**2)
    else:
        Llin = NS + (Lict-NS)*(CS-CE)/(CS-CICT)
        dLE_dCS = dLlin_dCS * (1.+c1+c2*Llin)
        dLE_dCICT = dLlin_dCICT * (1.+c1+c2*Llin)
        dLE_dCE = dLlin_dCE * (1.+c1+c2*Llin)
        dLE_dLict = dLlin_dLict * (1.+c1+c2*Llin)

        uncert = (dLE_dCS**2)*(av_noise**2) + \
                 (dLE_dCICT**2)*(ict_noise**2) + \
                 (dLE_dCE**2)*(noise**2) + \
                 (dLE_dLict**2)*(ict_random**2)
        
    return np.sqrt(uncert)

def get_sys(channel,uICT,Tict,CS,CE,CICT,NS,c1,c2,convT):
    """Get the systematic parts of the IR calibration uncertainty. Done per 
    scanline"""
    #
    # Gain part for all noise sources
    #
    Lict,uradTict = convT.t_to_rad_uncert(Tict,uICT)
    dLlin_dLict = (CS-CE)/(CS-CICT)

    #
    # If channel = 2,3 (11/12) then add non-linear part
    #
    if channel == 1:
        uncert = (dLlin_dLict**2)*(uradTict**2)
    else:
        Llin = NS + (Lict-NS)*(CS-CE)/(CS-CICT)
        dLE_dLict = dLlin_dLict * (1.+c1+c2*Llin)

        uncert = (dLE_dLict**2)*(uradTict**2)

    return np.sqrt(uncert)

def get_vars(ds,channel,convT,wlength,prt_threshold,ict_threshold,
             space_threshold,gac,cal,mask,out_prt=False,out_solza=False,
             out_time=False):
    """Get variables from xarray including smoothing and interpolation"""

    space = ds['space_counts'].values[:,channel]
    prt = ds["prt_counts"].values[:]
    ict = ds['ict_counts'].values[:,channel]
    ce = ds['counts'].values[:,:,channel-3]
    midpoint = ds['sun_zen'].shape[1]//2
    line_numbers = ds["scan_line_index"].data
    
    if out_solza:
        solza = ds['sun_zen'].values[:,midpoint]
    if out_time:
        time = (ds['times'].values[:] - 
                np.datetime64('1970-01-01T00:00:00Z'))/\
                np.timedelta64(1,'s')
        
    #
    # Set nan's to value to be caught by interpolation routines
    #
    gd = ~np.isfinite(prt)
    if np.sum(gd) > 0:
        prt[gd] = 0
    gd = ~np.isfinite(ict)
    if np.sum(gd) > 0:
        ict[gd] = 0
    gd = ~np.isfinite(space)
    if np.sum(gd) > 0:
        space[gd] = 0

    #
    # Removed old way of doing PRT indexing
    #
    # PRT index check
    #
    # PRTs. See reader.get_telemetry implementations.
    #
    #for offset in range(5):
    #    # According to the KLM Guide the fill value between PRT measurments is 0, but we search
    #    # for the first measurement gap using the threshold, because the fill value is in practice
    #    # not always exactly 0.
    #    if np.median(prt[(line_numbers - line_numbers[0]) % 5 == offset]) < prt_threshold:
    #        break
    #    else:
    #        raise IndexError("No PRT 0-index found!")
    #
    # get the PRT index, iprt equals to 0 corresponds to the measurement gaps
    # This can give wrong ICT temperatures
    #
    #iprt = (line_numbers - line_numbers[0] + 5 - offset) % 5
    
    # Get PRT mapping using new technique
    iprt = get_prt_nos(prt,prt_threshold,line_numbers,gac)
    
    #
    # Interpolate over bad prt values - from pygac calibrate_thermal
    #
    # fill measured values below threshold by interpolation
    #
    ifix = np.where(np.logical_and(iprt == 1, prt <= prt_threshold))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 1, prt > prt_threshold))
        if len(inofix[0]):
            prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])
        else:
            raise IndexError('No good prt1 data')

    ifix = np.where(np.logical_and(iprt == 2, prt <= prt_threshold))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 2, prt > prt_threshold))
        if len(inofix[0]):
            prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])
        else:
            raise IndexError('No good prt2 data')

    ifix = np.where(np.logical_and(iprt == 3, prt <= prt_threshold))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 3, prt > prt_threshold))
        if len(inofix[0]):
            prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])
        else:
            raise IndexError('No good prt3 data')

    ifix = np.where(np.logical_and(iprt == 4, prt <= prt_threshold))
    if len(ifix[0]):
        inofix = np.where(np.logical_and(iprt == 4, prt > prt_threshold))
        if len(inofix[0]):
            prt[ifix] = np.interp(ifix[0], inofix[0], prt[inofix])    
        else:
            raise IndexError('No good prt4 data')

    #
    # Convert to temperature
    #
    # calculate PRT temperature using equation (7.1.2.4-1) KLM Guide
    # Tprt = d0 + d1*Cprt + d2*Cprt^2 + d3*Cprt^3 + d4*Cprt^4
    # Note: First dimension of cal.d are the five coefficient indicees
    #
    tprt = np.polynomial.polynomial.polyval(prt, cal.d[:, iprt], tensor=False)

    #
    # Get interpolated values as done in pygac
    #
    tprt_interp = np.copy(tprt)
    zeros = iprt == 0
    nonzeros = np.logical_not(zeros)

    tprt_interp[zeros] = np.interp((zeros).nonzero()[0],
                            (nonzeros).nonzero()[0],
                            tprt[nonzeros])
    #
    # Interpolate over each PRT number
    #
    tprt1_interp = np.copy(tprt)
    zeros = (iprt == 0)|(iprt != 1)
    nonzeros = np.logical_not(zeros)

    tprt1_interp[zeros] = np.interp((zeros).nonzero()[0],
                            (nonzeros).nonzero()[0],
                            tprt[nonzeros])

    tprt2_interp = np.copy(tprt)
    zeros = (iprt == 0)|(iprt != 2)
    nonzeros = np.logical_not(zeros)

    tprt2_interp[zeros] = np.interp((zeros).nonzero()[0],
                            (nonzeros).nonzero()[0],
                            tprt[nonzeros])

    tprt3_interp = np.copy(tprt)
    zeros = (iprt == 0)|(iprt != 3)
    nonzeros = np.logical_not(zeros)

    tprt3_interp[zeros] = np.interp((zeros).nonzero()[0],
                            (nonzeros).nonzero()[0],
                            tprt[nonzeros])

    tprt4_interp = np.copy(tprt)
    zeros = (iprt == 0)|(iprt != 4)
    nonzeros = np.logical_not(zeros)

    tprt4_interp[zeros] = np.interp((zeros).nonzero()[0],
                            (nonzeros).nonzero()[0],
                            tprt[nonzeros])
    
    # Thresholds to flag missing/wrong data for interpolation
    # Remove masked data
    ict[mask] = 0
    space[mask] = 0
    zeros = ict < ict_threshold
    nonzeros = np.logical_not(zeros)
    try:
        ict[zeros] = np.interp((zeros).nonzero()[0],
                               (nonzeros).nonzero()[0],
                               ict[nonzeros])
    except ValueError: # 3b has no valid data
        raise IndexError('Channel 3b has no valid data')
    zeros = space < space_threshold
    nonzeros = np.logical_not(zeros)

    space[zeros] = np.interp((zeros).nonzero()[0],
                             (nonzeros).nonzero()[0],
                             space[nonzeros])
    
    #
    # Make averages and do using pygacs method at this point
    #
    weighting_function = np.ones(wlength, dtype=float) / wlength
    tprt_convolved = np.convolve(tprt_interp, weighting_function, "same")
    tprt1_convolved = np.convolve(tprt1_interp, weighting_function, "same")
    tprt2_convolved = np.convolve(tprt2_interp, weighting_function, "same")
    tprt3_convolved = np.convolve(tprt3_interp, weighting_function, "same")
    tprt4_convolved = np.convolve(tprt4_interp, weighting_function, "same")
    ict_convolved = np.convolve(ict, weighting_function, "same")
    space_convolved = np.convolve(space, weighting_function, "same")

    # take care of the beginning and end
    tprt_convolved[0:(wlength - 1) // 2] = tprt_convolved[(wlength - 1) // 2]
    tprt1_convolved[0:(wlength - 1) // 2] = tprt1_convolved[(wlength - 1) // 2]
    tprt2_convolved[0:(wlength - 1) // 2] = tprt2_convolved[(wlength - 1) // 2]
    tprt3_convolved[0:(wlength - 1) // 2] = tprt3_convolved[(wlength - 1) // 2]
    tprt4_convolved[0:(wlength - 1) // 2] = tprt4_convolved[(wlength - 1) // 2]
    ict_convolved[0:(wlength - 1) // 2] = ict_convolved[(wlength - 1) // 2]
    space_convolved[0:(wlength - 1) // 2] = space_convolved[(wlength - 1) // 2]
    tprt_convolved[-(wlength - 1) // 2:] = tprt_convolved[-((wlength + 1) // 2)]
    tprt1_convolved[-(wlength - 1) // 2:] = tprt1_convolved[-((wlength + 1) // 2)]
    tprt2_convolved[-(wlength - 1) // 2:] = tprt2_convolved[-((wlength + 1) // 2)]
    tprt3_convolved[-(wlength - 1) // 2:] = tprt3_convolved[-((wlength + 1) // 2)]
    tprt4_convolved[-(wlength - 1) // 2:] = tprt4_convolved[-((wlength + 1) // 2)]
    ict_convolved[-(wlength - 1) // 2:] = ict_convolved[-((wlength + 1) // 2)]
    space_convolved[-(wlength - 1) // 2:] = space_convolved[-((wlength + 1) // 2)]

    if out_prt:
        if out_solza:
            if out_time:
                return space_convolved,ict_convolved,ce,tprt_convolved,\
                    tprt1_convolved,tprt2_convolved,tprt3_convolved,\
                    tprt4_convolved,solza,time
            else:
                return space_convolved,ict_convolved,ce,tprt_convolved,\
                    tprt1_convolved,tprt2_convolved,tprt3_convolved,\
                    tprt4_convolved,solza
        else:
            return space_convolved,ict_convolved,ce,tprt_convolved,\
                tprt1_convolved,tprt2_convolved,tprt3_convolved,\
                tprt4_convolved
    else:
        if out_solza:
            if out_time:
                return space_convolved,ict_convolved,ce,tprt_convolved,solza,\
                time
            else:
                return space_convolved,ict_convolved,ce,tprt_convolved,solza
        else:
            return space_convolved,ict_convolved,ce,tprt_convolved

def get_gainval(time,avhrr,ict1,ict2,ict3,ict4,CS,CICT,CE,NS,bad_scan,convT,
                window,calculate=False):
    """Estimate gain value at smallest stdev point in orbit either from
    file or estimate it from data"""
    #
    # Open file containing 3.7mu gain value and interpolate over time
    #
    # Now stored on zenodo
    #
    #coef_file = files("pygac") / "data/{0}_uncert.nc".format(avhrr)
    coef_file = "https://zenodo.org/records/15482385/files/{0}_uncert.nc#mode=bytes".format(avhrr)
    try:
        d = xr.open_dataset(coef_file)
    except FileNotFoundError:
        raise Exception("ERROR: Gain can not be determined because zenodo not available")
    intime = d["time_gain"].values[:]
    ingain = d["gain"].values[:]
    d.close()
    
    if time > intime[-1] and not calculate:
        raise Exception('ERROR: file time > last time with max gain values (HRPR/LAC)')
        
    timediff = (intime-time)/np.timedelta64(1,'s')
    timediff = np.abs(timediff)
    timediff_min = timediff.min()
    #
    # Within a day at worst or if there is no 3B data available (3A/3B switch
    # case)
    #
    # Check Earth counts in window for 3.7 micron as if nan then in S3/3B
    #
    
    if not calculate:
        #
        # HRPT/LAC
        #
        pos = np.nonzero(timediff == timediff_min)[0][0]
        return ingain[pos]        
    elif timediff_min < 86400.:
        #
        # Take precalc value from file
        #
        pos = np.nonzero(timediff == timediff_min)[0][0]
        return ingain[pos]        
    else:
        #
        # No nearby gain estimate or force calculate so calculate
        # from data using min std of prts
        #
        gd = (bad_scan == 0)
        prt1 = prt1[gd]
        prt2 = prt2[gd]
        prt3 = prt3[gd]
        prt4 = prt4[gd]
        Lict = Lict[gd]
        CS = CS[gd]
        CICT = CICT[gd]
        #
        # Get stdev
        #
        X = np.zeros((len(prt1),4))
        X[:,0] = prt1
        X[:,1] = prt2
        X[:,2] = prt3
        X[:,3] = prt4
        stdev = np.std(X,axis=1)
        pos = np.nonzero(stdev == stdev.min())[0][0]
        Lict = convT.t_to_rad(np.mean(X,axis=1))
        gain = (Lict-NS)/(CS-CICT)
        return gain[pos]

def get_pixel(Lict,CS,CE,CICT,NS,c0,c1,c2):
    """Get radiance of pixel using calibration"""

    Llin = NS + (Lict-NS)*(CS-CE)/(CS-CICT)
    if NS != 0.:
        LE = Llin + c0 + c1*Llin + c2*Llin*Llin
        return LE
    else:
        return Llin

def get_uncert_parameter_thresholds(vischans=False):
    """Return required constants for IR and visible uncertainty cases. Single
    point so any changes will be correctly applied across both"""
    #
    # Define averaging kernel based on value in noaa.py
    # Also set PRT uncertainty components and thresholds
    #
    if vischans:
        window = 51
        solar_contam_threshold = 0.05
        solar_contam_sza_threshold = 102.
        return window,solar_contam_threshold,solar_contam_sza_threshold
    else:
        window = 51
        prt_bias = 0.01
        prt_sys = 0.1
        prt_threshold = 50
        ict_threshold = 100
        space_threshold = 100

        return window, prt_bias, prt_sys, prt_threshold, ict_threshold, \
            space_threshold

def get_solar_from_file(avhrr_name,ds):
    """Read in possible solar contamination times from uncertainty files.
    Only used for LAC/HRPT data where there is not enough data to detect
    possible solar contamination so GAC estimates are used"""

    #
    # Get times in seconds from
    #
    time = (ds['times'].values[:] - \
            np.datetime64('1970-01-01T00:00:00Z'))/\
            np.timedelta64(1,'s')
    
    #
    # Read file
    #
    # Now stored on zenodo
    #
    # coef_file = files("pygac") / "data/{0}_uncert.nc".format(avhrr_name)
    coef_file = "https://zenodo.org/records/15482385/files/{0}_uncert.nc#mode=bytes".format(avhrr_name)
    try:
        d = xr.open_dataset(coef_file,decode_times=False)
    except FileNotFoundError:
        raise Exception("ERROR: Solar data can mot be determined because zenodo not available")
    solar_start_time_1 = d["gain1_solar_start"].values[:]
    solar_stop_time_1 = d["gain1_solar_stop"].values[:]
    solar_start_time_2 = d["gain2_solar_start"].values[:]
    solar_stop_time_2 = d["gain2_solar_stop"].values[:]
    d.close()
    
    gd = np.isfinite(solar_start_time_1)&np.isfinite(solar_stop_time_1)
    solar_start_time_1 = solar_start_time_1[gd]
    solar_stop_time_1 = solar_stop_time_1[gd]
    gd = np.isfinite(solar_start_time_2)&np.isfinite(solar_stop_time_2)
    solar_start_time_2 = solar_start_time_2[gd]
    solar_stop_time_2 = solar_stop_time_2[gd]
    #
    # Match to times in file
    #
    min_solar_1 = -1
    max_solar_1 = -1
    for i in range(len(solar_start_time_1)):
        gd = (time >= solar_start_time_1[i])&(time <= solar_stop_time_1[i])
        if np.sum(gd) > 0:
            index = np.arange(len(time)).astype(dtype=np.int32)
            index = index[gd]
            min_solar_1 = index[0]
            max_solar_1 = index[-1]
    min_solar_2 = -1
    max_solar_2 = -1
    for i in range(len(solar_start_time_2)):
        gd = (time >= solar_start_time_2[i])&(time <= solar_stop_time_2[i])
        if np.sum(gd) > 0:
            index = np.arange(len(time)).astype(dtype=np.int32)
            index = index[gd]
            min_solar_2 = index[0]
            max_solar_2 = index[-1]

    #
    # Return solar location if available
    #
    return min_solar_1, max_solar_1, min_solar_2, max_solar_2

def ir_uncertainty(ds,mask,plot=False,plotmax=None):
    """Create the uncertainty components for the IR channels. These include
    
    1) Random 
          a) Noise
          b) Digitisation
          c) ICT PRT Noise
    2) Systematic
          a) ICT Temperature uncertainty
          b) PRT Bias
          c) Calibration coefs/measurement equation uncertainty
    
    Inputs:
          ds : Input xarray dataset containing data for calibration
        mask : pygac mask from reader
    Outputs:
      uncert : xarray dataset containing random and systematic uncertainty
               components
    """

    #
    # Get parameters for kernals/uncertainties
    #
    window, prt_bias, prt_sys, prt_threshold, ict_threshold, \
        space_threshold = get_uncert_parameter_thresholds()
    
    if ds['channels'].values.shape[1] == 409:
        gacdata = True
        scale = 1
    else:
        gacdata = False
        scale = 3

    avhrr_name = ds.attrs["spacecraft_name"]

    #
    # Get calibration coefficients
    #
    cal = Calibrator(
        ds.attrs["spacecraft_name"])
    NS_1 = cal.space_radiance[0]
    NS_2 = cal.space_radiance[1]
    NS_3 = cal.space_radiance[2]
    c0_1 = cal.b[0,0]
    c1_1 = cal.b[0,1]
    c2_1 = cal.b[0,2]
    c0_2 = cal.b[1,0]
    c1_2 = cal.b[1,1]
    c2_2 = cal.b[1,2]
    c0_3 = cal.b[2,0]
    c1_3 = cal.b[2,1]
    c2_3 = cal.b[2,2]

    # Is the twelve micron channel there
    if ds.attrs["spacecraft_name"] == "tirosn" or \
       ds.attrs["spacecraft_name"] == "noaa06" or \
       ds.attrs["spacecraft_name"] == "noaa08" or \
       ds.attrs["spacecraft_name"] == "noaa10":
        twelve_micron = False
    else:
        twelve_micron = True

    # Temperature to radiance etc.
    convT1 = convBT(cal,0)
    convT2 = convBT(cal,1)
    if twelve_micron:
        convT3 = convBT(cal,2)

    #
    # Get variables for 10 sampled case
    #
    total_space = ds['total_space_counts'].values[:,:,:]
    total_ict = ds['total_ict_counts'].values[:,:,:]
    #
    # Noise elements
    #
    noise1,noise2,noise3,av_noise1,av_noise2,av_noise3,\
        av_ict_noise1,av_ict_noise2,av_ict_noise3,bad_scan \
        = get_noise(total_space,total_ict,window,twelve_micron)

    #
    # Get variables used on the calibration
    #
    nfigure=1
    if plot:
        import matplotlib.pyplot as plt
        plt.figure(nfigure)
    else:
        plt = None

    CS_1,CICT_1,CE_1,Tict,ict1,ict2,ict3,ict4 = get_vars(ds,0,convT1,
                                                         window,
                                                         prt_threshold,
                                                         ict_threshold,
                                                         space_threshold,
                                                         gacdata,
                                                         cal,
                                                         mask,
                                                         out_prt=True)
    if plot:
        plt.subplot(131)
        plt.plot(np.arange(len(CS_1)),CS_1,',')
        plt.title('3.7$\mu$m')
        plt.ylabel('Space Cnts')
        plt.xlabel('Scanline')

    CS_2,CICT_2,CE_2,Tict = get_vars(ds,1,convT2,window,prt_threshold,
                                     ict_threshold,
                                     space_threshold,
                                     gacdata,cal,mask)
    if plot:
        plt.subplot(132)
        plt.plot(np.arange(len(CS_2)),CS_2,',')
        plt.title('11$\mu$m')
        plt.ylabel('Space Cnts')
        plt.xlabel('Scanline')
        
    if twelve_micron:
        CS_3,CICT_3,CE_3,Tict = get_vars(ds,2,convT3,window,
                                         prt_threshold,
                                         ict_threshold,
                                         space_threshold,
                                         gacdata,cal,mask)
        if plot:
            plt.subplot(133)
            plt.plot(np.arange(len(CS_3)),CS_3,',')
            plt.title('12$\mu$m')
            plt.ylabel('Space Cnts')
            plt.xlabel('Scanline')
    if plot:
        plt.tight_layout()

    solar_flag = np.zeros(CE_2.shape[0],dtype=np.uint8)
    if gacdata:
        #
        # See if solar contamination present
        # Only for GAC data
        # Possible at 2 locations
        #
        min_solar_1, max_solar_1, peak_solar_1, solar_solza_1,\
        min_solar_2, max_solar_2, peak_solar_2, solar_solza_2 = \
            find_solar(ds,mask,convT1,plot=plot,out_time=False)
        if min_solar_1 is None and max_solar_1 is None:
            min_solar_1 = -1
            max_solar_1 = -1
        if min_solar_2 is None and max_solar_2 is None:
            min_solar_2 = -1
            max_solar_2 = -1
    else:
        #
        # Find if stored solar contamination is present
        #
        min_solar_1, max_solar_1, min_solar_2, max_solar_2 = \
            get_solar_from_file(avhrr_name,ds)
    #
    # Set solar flag
    #
    if min_solar_1 >= 0 and max_solar_1 >= 0:
        solar_flag[min_solar_1:max_solar_1+1] = 1
    if min_solar_2 >= 0 and max_solar_2 >= 0:
        solar_flag[min_solar_2:max_solar_2+1] = 1

    #
    # Systematic components - uICT from gain variation in 3.7mu channel
    #
    gd = np.isfinite(ds['times'].values)
    time = ds['times'].values[gd][0]
    #
    # Only redo calculation if gac data
    #
    gain_37 = get_gainval(time,avhrr_name,ict1,ict2,ict3,ict4,CS_1,CICT_1,
                          CE_1,0.,bad_scan,convT1,window,calculate=gacdata)
    uICT,nfigure = get_uICT(gain_37,CS_1,CICT_1,Tict,0.,convT1,bad_scan,
                            solar_flag,window,plt=plt,nfigure=nfigure)
    #
    # Loop round scanlines
    #
    bt_rand_37 = np.zeros(CE_2.shape,dtype=CE_2.dtype)
    bt_rand_11 = np.zeros(CE_2.shape,dtype=CE_2.dtype)
    bt_rand_12 = np.zeros(CE_2.shape,dtype=CE_2.dtype)
    if not twelve_micron:
        bt_rand_12[:,:] = np.nan
    bt_sys_37 = np.zeros(CE_2.shape,dtype=CE_2.dtype)
    bt_sys_11 = np.zeros(CE_2.shape,dtype=CE_2.dtype)
    bt_sys_12 = np.zeros(CE_2.shape,dtype=CE_2.dtype)
    if not twelve_micron:
        bt_sys_12[:,:] = np.nan
    urand_37 = np.zeros(CE_2.shape,dtype=CE_2.dtype)
    urand_11 = np.zeros(CE_2.shape,dtype=CE_2.dtype)
    urand_12 = np.zeros(CE_2.shape,dtype=CE_2.dtype)
    if not twelve_micron:
        urand_12[:,:] = np.nan
    uratio_37 = np.zeros(CE_2.shape,dtype=CE_2.dtype)
    uratio_11 = np.zeros(CE_2.shape,dtype=CE_2.dtype)
    uratio_12 = np.zeros(CE_2.shape,dtype=CE_2.dtype)
    if not twelve_micron:
        uratio_12[:,:] = np.nan
    for i in range(len(CS_2)):
        #
        # Check for bad scanlines
        #
        if bad_scan[i] == 1:
            bt_rand_37[i,:] = np.nan
            bt_rand_11[i,:] = np.nan
            if twelve_micron:
                bt_rand_12[i,:] = np.nan
            bt_sys_37[i,:] = np.nan
            bt_sys_11[i,:] = np.nan
            if twelve_micron:
                bt_sys_12[i,:] = np.nan
            continue

        #
        # Get uncertainty in ICT temperature from PRT measurements
        #
        ict_random1, ict_sys1, Lict_1 = \
            get_ict_uncert(Tict[i],prt_bias,prt_sys,uICT[i],convT1)
        ict_random2, ict_sys2, Lict_2 = \
            get_ict_uncert(Tict[i],prt_bias,prt_sys,uICT[i],convT2)
        if twelve_micron:
            ict_random3, ict_sys3, Lict_3 = \
                get_ict_uncert(Tict[i],prt_bias,prt_sys,uICT[i],convT3)

        #
        # get pixel radiance
        #
        rad_37 = get_pixel(Lict_1,CS_1[i],
                           CE_1[i,:],CICT_1[i],0.,0.,0.,0.)
        rad_11 = get_pixel(Lict_2,CS_2[i],
                           CE_2[i,:],CICT_2[i],NS_2,c0_2,c1_2,c2_2)
        if twelve_micron:
            rad_12 = get_pixel(Lict_3,CS_3[i],
                               CE_3[i,:],CICT_3[i],NS_3,c0_3,c1_3,c2_3)

        #
        # Get noise in radiance space
        #
        rad_noise_37 = get_random(1,noise1,av_noise1,av_ict_noise1,
                                  ict_random1,Lict_1,CS_1[i],
                                  CE_1[i,:],CICT_1[i],0.,0.,0.)
        rad_noise_11 = get_random(2,noise2,av_noise2,av_ict_noise2,
                                  ict_random2,Lict_2,CS_2[i],
                                  CE_2[i,:],CICT_2[i],NS_2,c1_2,c2_2)
        if twelve_micron:
            rad_noise_12 = get_random(3,noise3,av_noise3,
                                      av_ict_noise3,
                                      ict_random3,Lict_3,CS_3[i],
                                      CE_3[i,:],CICT_3[i],NS_3,c1_3,c2_3)

        #
        # Convert to BT space uncertainty
        #
        T,bt_rand_37[i,:] = convT1.rad_to_t_uncert(rad_37,rad_noise_37)
        T,bt_rand_11[i,:] = convT2.rad_to_t_uncert(rad_11,rad_noise_11)
        if twelve_micron:
            T,bt_rand_12[i,:] = convT3.rad_to_t_uncert(rad_12,rad_noise_12)

        #
        # Get systematic uncertainty through the measurement equation
        # Note measurement equation uncertainty set to 0.5K@300
        #
        rad_sys_37 = get_sys(1,ict_sys1,Tict[i],CS_1[i],
                             CE_1[i,:],CICT_1[i],0.,0.,0.,convT1)
        rad_sys_11 = get_sys(2,ict_sys2,Tict[i],CS_2[i],
                             CE_2[i,:],CICT_2[i],NS_2,c1_2,c2_2,
                             convT2)
        if twelve_micron:
            rad_sys_12 = get_sys(3,ict_sys3,Tict[i],CS_3[i],
                                 CE_3[i,:],CICT_3[i],NS_3,c1_3,
                                 c2_3,convT3)

        T,bt_sys_37[i,:] = convT1.rad_to_t_uncert(rad_37,rad_sys_37)
        T,bt_sys_11[i,:] = convT2.rad_to_t_uncert(rad_11,rad_sys_11)
        if twelve_micron:
            T,bt_sys_12[i,:] = convT3.rad_to_t_uncert(rad_12,rad_sys_12)
            
        #
        # Add 0.5/sqrt(3.)K for measurement equation uncertainty
        # Also get ratio of ICT/total uncertainty
        # For 3.7mu channel add radiance due to Planck function
        # 0.5K @ 300K
        rad1 = convT1.t_to_rad(300.-0.5/np.sqrt(3.))
        rad2 = convT1.t_to_rad(300.+0.5/np.sqrt(3.))
        delta_rad = (rad2-rad1)/2.
        T,new_37_uncert = convT1.rad_to_t_uncert(rad_37,delta_rad)
        #tot_sys_37 = np.sqrt(bt_sys_37[i,:]**2+0.5**2/3.)
        tot_sys_37 = np.sqrt(bt_sys_37[i,:]**2+new_37_uncert**2)
        #
        # Just add 0.5K (11/12)
        #
        tot_sys_11 = np.sqrt(bt_sys_11[i,:]**2+0.5**2/3.)
        if twelve_micron:
            tot_sys_12 = np.sqrt(bt_sys_12[i,:]**2+0.5**2/3.)
        uratio_37[i,:] = bt_sys_37[i,:] / tot_sys_37 
        uratio_11[i,:] = bt_sys_11[i,:] / tot_sys_11 
        if twelve_micron:
            uratio_12[i,:] = bt_sys_12[i,:] / tot_sys_12
        bt_sys_37[i,:] = tot_sys_37
        bt_sys_11[i,:] = tot_sys_11
        if twelve_micron:
            bt_sys_12[i,:] = tot_sys_12

    if plot:
        if twelve_micron:
            nfigure = nfigure+1
            plt.figure(nfigure)
            plt.subplot(231)
            plt.hist(bt_rand_37.flatten(),bins=100)
            plt.title('3.7$\mu$m')
            
            plt.subplot(232)
            plt.hist(bt_rand_11.flatten(),bins=100)
            plt.title('11$\mu$m (Random)')

            plt.subplot(233)
            plt.hist(bt_rand_12.flatten(),bins=100)
            plt.title('12$\mu$m')
            
            plt.subplot(234)
            plt.hist(bt_sys_37.flatten(),bins=100)
            plt.title('3.7$\mu$m')
            
            plt.subplot(235)
            plt.hist(bt_sys_11.flatten(),bins=100)
            plt.title('11$\mu$m (Systematic)')
            plt.xlabel('Uncertainty / K')

            plt.subplot(236)
            plt.hist(bt_sys_12.flatten(),bins=100)
            plt.title('12$\mu$m')
            plt.tight_layout()

            nfigure = nfigure+1
            plt.figure(nfigure)
            plt.subplot(231)
            if plotmax is None:
                im=plt.imshow(bt_rand_37)
            else:
                gd = (bt_rand_37 > plotmax[0])
                image = np.copy(bt_rand_37)
                image[gd] = np.nan
                im=plt.imshow(image)
            plt.colorbar(im)
            plt.title('3.7$\mu$m')

            plt.subplot(232)
            if plotmax is None:
                im=plt.imshow(bt_rand_11)
            else:
                gd = (bt_rand_11 > plotmax[1])
                image = np.copy(bt_rand_11)
                image[gd] = np.nan
                im=plt.imshow(image)
            plt.colorbar(im)
            plt.title('11$\mu$m (Random)')

            plt.subplot(233)
            if plotmax is None:
                im=plt.imshow(bt_rand_12)
            else:
                gd = (bt_rand_12 > plotmax[2])
                image = np.copy(bt_rand_12)
                image[gd] = np.nan
                im=plt.imshow(image)
            plt.colorbar(im)
            plt.title('12$\mu$m')

            plt.subplot(234)
            if plotmax is None:
                im=plt.imshow(bt_sys_37)
            else:
                gd = (bt_sys_37 > plotmax[3])
                image = np.copy(bt_sys_37)
                image[gd] = np.nan
                im=plt.imshow(image)
            plt.colorbar(im)
            plt.title('3.7$\mu$m')

            plt.subplot(235)
            if plotmax is None:
                im=plt.imshow(bt_sys_11)
            else:
                gd = (bt_sys_11 > plotmax[4])
                image = np.copy(bt_sys_11)
                image[gd] = np.nan
                im=plt.imshow(image)
            plt.colorbar(im)
            plt.title('11$\mu$m (Systematic)')

            plt.subplot(236)
            if plotmax is None:
                im=plt.imshow(bt_sys_12)
            else:
                gd = (bt_sys_12 > plotmax[5])
                image = np.copy(bt_sys_12)
                image[gd] = np.nan
                im=plt.imshow(image)
            plt.colorbar(im)
            plt.title('12$\mu$m')
            plt.tight_layout()

            nfigure = nfigure+1
            plt.figure(nfigure)
            plt.subplot(131)
            im=plt.imshow(uratio_37)
            plt.colorbar(im)
            plt.title('3.7$\mu$m')
            plt.subplot(132)
            im=plt.imshow(uratio_11)
            plt.colorbar(im)
            plt.title('11$\mu$m')
            plt.subplot(133)
            im=plt.imshow(uratio_12)
            plt.colorbar(im)
            plt.title('12$\mu$m')
            plt.tight_layout()
        else:
            nfigure = nfigure+1
            plt.figure(nfigure)
            plt.subplot(221)
            plt.hist(bt_rand_37.flatten(),bins=100)
            plt.title('3.7$\mu$m (Random)')
            
            plt.subplot(222)
            plt.hist(bt_rand_11.flatten(),bins=100)
            plt.title('11$\mu$m (Random)')

            plt.subplot(223)
            plt.hist(bt_sys_37.flatten(),bins=100)
            plt.title('3.7$\mu$m (Systematic)')
            plt.xlabel('Uncertainty / K')
            
            plt.subplot(224)
            plt.hist(bt_sys_11.flatten(),bins=100)
            plt.title('11$\mu$m (Systematic)')
            plt.xlabel('Uncertainty / K')

            plt.tight_layout()
        plt.show()
    #
    # Output uncertainties
    #
    random = np.zeros((bt_rand_11.shape[0],bt_rand_11.shape[1],3))
    systematic = np.zeros((bt_rand_11.shape[0],bt_rand_11.shape[1],3))
    uratio = np.zeros((bt_rand_11.shape[0],bt_rand_11.shape[1],3),
                      dtype=np.uint8)
    uflags = np.zeros((bt_rand_11.shape[0]),dtype=np.uint8)

    random[:,:,0] = bt_rand_37
    random[:,:,1] = bt_rand_11
    if twelve_micron:
        random[:,:,2] = bt_rand_12
    else:
        random[:,:,2] = np.nan
    systematic[:,:,0] = bt_sys_37
    systematic[:,:,1] = bt_sys_11
    if twelve_micron:
        systematic[:,:,2] = bt_sys_12
    else:
        systematic[:,:,2] = np.nan
    #
    # Ratio for channel-to-channel covariance as ubyte
    #
    gd = (uratio_37 < 0.)
    uratio_37[gd] = 0.
    gd = (uratio_37 > 1.)
    uratio_37[gd] = 1.
    gd = (uratio_11 < 0.)
    uratio_11[gd] = 0.
    gd = (uratio_11 > 1.)
    uratio_11[gd] = 1.
    gd = np.isfinite(uratio_37)
    uratio[gd,0] = (uratio_37[gd]*255).astype(dtype=np.uint8)
    gd = np.isfinite(uratio_11)
    uratio[gd,1] = (uratio_11[gd]*255).astype(dtype=np.uint8)
    if twelve_micron:
        gd = (uratio_12 < 0.)
        uratio_12[gd] = 0.
        gd = (uratio_12 > 1.)
        uratio_12[gd] = 1.
        gd = np.isfinite(uratio_12)
        uratio[gd,2] = (uratio_12[gd]*255).astype(dtype=np.uint8)
    else:
        uratio[:,:,2] = 0

    #
    # Flags
    #
    gd = (bad_scan == 1)
    uflags[gd] = 1
    gd = (solar_flag == 1)
    uflags[gd] = (uflags[gd]|2)

    time = (ds["times"].values - np.datetime64("1970-01-01 00:00:00"))/\
           np.timedelta64(1,'s')
    time_da = xr.DataArray(time,dims=["times"],attrs={"long_name":"scanline time",
                                                     "units":"seconds since 1970-01-01"})
    across_da = xr.DataArray(np.arange(random.shape[1]),dims=["across_track"])
    ir_channels_da = xr.DataArray(np.array([3,4,5]),dims=["ir_channels"])
    random_da = xr.DataArray(random,dims=["times","across_track","ir_channels"],
                             attrs={"long_name":"Random uncertainties","units":"K"})
    sys_da = xr.DataArray(systematic,dims=["times","across_track","ir_channels"],
                          attrs={"long_name":"Systematic uncertainties","units":"K"})

    uratio_da = xr.DataArray(uratio,dims=["times","across_track","ir_channels"],
                             attrs={"long_name":"Channel-to-channel covariance  ratio",
                                    "_FillValue":0})

    uflags_da = xr.DataArray(uflags,dims=["times"],
                             attrs={"long_name":"Uncertainty flags (bit 1==bad space view (value=1), bit 2==solar contamination (value=2)"})

    uncertainties = xr.Dataset(dict(times=time_da,across_track=across_da,
                                    ir_channels=ir_channels_da,
                                    random=random_da,systematic=sys_da,
                                    chan_covar_ratio=uratio_da,
                                    uncert_flags=uflags_da))

    return uncertainties

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('--plot',action='store_true')
    parser.add_argument('--plotmax',type=float,nargs=6)
    parser.add_argument('--solar_contam',action='store_true')
    parser.add_argument('--oname')
    parser.add_argument('--write_test',nargs=2)
    
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
    if args.solar_contam:
        if args.oname is None:
            raise Exception('ERROR: if solar data wanted must use --oname')
        if args.plot:
            import matplotlib.pyplot as plt
        tmin_1,tmax_1,tpeak_1,out_solza_1,\
        tmin_2,tmax_2,tpeak_2,out_solza_2,\
        gtime,min_gain = \
            find_solar(ds,mask,convT=None,outgain=True,plot=args.plot,
                       out_time=True)
        with open(args.oname,'a') as fp:
            if tmin_1 is not None and tmin_2 is not None:
                fp.write('{0:16.9e} {1:16.9e} {2:16.9e} {3:9.4f} {4:16.9e} {5:16.9e} {6:16.9e} {7:9.4f} {8:16.9e} {9:9.7e}\n'.
                     format(tmin_1,tmax_1,tpeak_1,out_solza_1,
                            tmin_2,tmax_2,tpeak_2,out_solza_2,
                            gtime,min_gain))
            elif tmin_1 is not None and tmin_2 is None:
                fp.write('{0:16.9e} {1:16.9e} {2:16.9e} {3:9.4f} {4:16.9e} {5:16.9e} {6:16.9e} {7:9.4f} {8:16.9e} {9:9.7e}\n'.
                     format(tmin_1,tmax_1,tpeak_1,out_solza_1,
                            -1.,-1.,-1.,-1.,
                            gtime,min_gain))
            elif tmin_1 is None and tmin_2 is not None:
                fp.write('{0:16.9e} {1:16.9e} {2:16.9e} {3:9.4f} {4:16.9e} {5:16.9e} {6:16.9e} {7:9.4f} {8:16.9e} {9:9.7e}\n'.
                     format(-1.,-1.,-1.,-1.,
                            tmin_2,tmax_2,tpeak_2,out_solza_2,
                            gtime,min_gain))
            elif tmin_1 is None and tmin_2 is None:
                fp.write('{0:16.9e} {1:16.9e} {2:16.9e} {3:9.4f} {4:16.9e} {5:16.9e} {6:16.9e} {7:9.4f} {8:16.9e} {9:9.7e}\n'.
                     format(-1.,-1.,-1.,-1.,
                            -1.,-1.,-1.,-1.,
                            gtime,min_gain))
    elif args.write_test is not None:
        #
        # Write data for testing of uncertainty code
        #
        if ds.attrs['midnight_scanline'] is None:
            ds.attrs['midnight_scanline'] = -1
        ds.to_netcdf(args.write_test[0])
        X = np.zeros((len(mask)),dtype=np.int8)
        X[mask] = 1
        np.savetxt(args.write_test[1],X,fmt="%d")
    else:
        uncert = ir_uncertainty(ds,mask,plot=args.plot,plotmax=args.plotmax)
        print(uncert)

