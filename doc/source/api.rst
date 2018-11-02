The :mod:`pygac` API
====================

*pygac* interface consists of a number of python modules. The following schematic shows a general structure and the processing flow of *pygac*.   

It must be noted that *pygac* expects Level 1b file to contain normal GAC header and data records, the format of which are mentioned in the official NOAA POD and KLM Data User Guides. The user should not prepend any other header (e.g. when downloading GAC data from CLASS archive etc) to the L1b file. In the first pre-processing step, *pygac* determines whether the GAC data comes from the second (i.e. NOAA-14 and before) or the third generation (NOAA-15 and onwards) AVHRR instrument by "pygac-run".
This is done by reading the first three bytes of the data set. If they contain the any of the following values, ["CMS", "NSS", "UKM", "DSS"], then the KLM reader from "gac_klm.py" file is invoked, otherwise the POD reader is invoked (gac_pod.py).


GAC POD reader
--------------

.. automodule:: pygac.gac_pod
   :members:
   :undoc-members:
   
As the format of GAC POD header is changed twice (once in 1992 and again in 1994), there are currently three readers integrated in the *pygac* to read POD header.

In case of POD family of satellites, the first data record often doesn't start with first scan line number. In fact, the latter often has unambiguous value and may not be the continuation of the last record number (as expected due to overlap).

Here is an example showing scanlines numbers for the first and the last 3 lines from one random NOAA-11 orbit:

[12804 3 4 ..., 12806 12807 12808]
You can see that the first two scanlines are missing and the orbit is actually starting with scanline number 3.

*pygac* would rearrange the entire orbit in the increasing order of the scanline numbers.

Sometimes the scanlines have erroneous time information. Here is an example of time stamp from few consecutive scanlines from a radom orbit from NOAA-11::

       [datetime.datetime(2071, 4, 28, 15, 44, 46, 329000)],
       [datetime.datetime(2067, 4, 27, 15, 44, 46, 336000)],
       [datetime.datetime(2065, 4, 25, 15, 44, 46, 328000)],
       [datetime.datetime(2059, 4, 22, 15, 44, 46, 331000)],
       [datetime.datetime(2057, 4, 22, 15, 27, 17, 757000)],
       [datetime.datetime(2055, 4, 20, 15, 27, 17, 753000)],
       [datetime.datetime(2069, 7, 20, 14, 17, 19, 377000)],
       [datetime.datetime(1981, 7, 25, 14, 17, 12, 379000)],
       [datetime.datetime(1987, 7, 27, 14, 17, 12, 380000)],
       [datetime.datetime(1991, 7, 29, 14, 17, 11, 359000)],
       [datetime.datetime(1993, 7, 29, 14, 17, 11, 360000)],

It is evident that year, month values are jumping and have non-sense values.

*pygac* would handle this in either of the following two ways. If the first scanline has valid time stamp, *pygac* would use this to calculate time stamp for remaining scanlines assuming the scanning rate of 0.5 sec/scanline and using actual scanline numbers (to account for any potential gaps in the data). If the first scanline doesn't have valid time stamp or missing altogether, then *pygac* computes the median value of valid time stamps in the data and extrapolates it based on scanline number and scanning rate to compute new time stamps for all scanlines.

Whenever possible, *pygac* uses RPY corrections along with other orbital parameters to compute accurate satellite location (e.g. instead of assuming constant altitude). However, RPY corrections are not available for all NOAA satellites. In case of the majority of the POD family satellites, these corrections are set to zero.


GAC KLM reader
--------------

.. automodule:: pygac.gac_klm
   :members:
   :undoc-members:


Computation of geolocation
--------------

.. automodule:: pygac.gac_klm
   :members:
   :undoc-members:


Each GAC row has total 409 pixels. But lat-lon values are provided for every eigth pixel starting from pixel 5 and ending at pixel 405. Using Numpy, Scipy and Pyresample packages, inter- and extrapolation is carried out to obtain geolocation for each pixel along the scanline.

If the GAC data belongs to POD family, then clock drift errors are used to adjust existing Lat-Lon information. Here, *pygac* makes use of PyOrbital package, which is a part of PyTroll suite of Python interface developed to process meteorological satellite data (further information here: http://www.pytroll.org/ and https://github.com/mraspaud/pyorbital). *pygac* interpolates the clock offset and adjusts the nominal scan times to the actual scan times. Since the geolocation was computed using the nominal scan times, *pygac* interpolates the latitudes and longitudes to the actual scan times using spherical linear interpolation, aka slerp. However, in the case of a clock drift error greater than the scan rate of the dataset, the latitude and longitude for each pixel of the scan lines that cannot have an interpolated geolocation (typically at the start or end of the dataset) are recomputed. This is done using pyorbital, which in turn uses TLEs to compute the position of the satellite at each scan time and the instrument geometry compute the longitude and latitude of each pixel of the dataset. Since this operation can be quite costly, the interpolation is prefered whenever possible.


GAC calibration/inter-calibration
--------------

.. automodule:: pygac.gac_klm
   :members:
   :undoc-members:

*pygac* currently supports calibration of all GAC data from AVHRRs onboard
 NOAA-7 and onwards, including Metop satellites.

At present, updated calibration coefficients provided by Andrew Heidinger (NOAA) under SCOPE-CM project are applied for all satellites. The solar channel calibration (Channels 1 and 2, and Channel 3a if available) takes into account inter-satellite differences and is derived using amalgamation of different calibration references including the most recent MODIS Collection 6 data, in-situ targets, and simultaneous nadir observations. The detailed methodology is presented in Heidinger et al. (2010). The resulting (inter)calibration coefficients are of highest quality and follow the Global Climate Observing System (GCOS) standards for deriving fundamental climate data records.

The reflectances are normalized by a factor (as a function of day of a year) to account for changing Earth-Sun distance. However, it is left to the potential user to apply further normalization using cosine of solar zenith angle (depending on application in question).


The thermal channel intercalibration is done from scratch, starting from obtaining Platinum Resistance Thermometer (PRT), space and Internal Calibration Target (ICT, blackbody) counts. ICT temperatures are obtained from PRT counts based on coefficients provided in the POD and KLM Data User Guides (Kidwell, 2000). For each thermal channel, a smoothing window of 51 successive PRT, ICT and space counts is used to obtain robust gain values and to dampen undue high frequency fluctuations in the count data (Trishchenko, 2002). This window size is easily configurable to suit user needs. Section 7.1.2.5 of KLM Data User Guide presents the summary of equations implemented in *pygac* to calibrate thermal channels (http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/c7/sec7-1.htm), including non-linearity correction (Walton et al. 1998).

The PRT readings are supposed to be present in a specific order, for example, the reset value followed by four readings from the PRTs. However, this may not always be the case for orbits that contain data gaps or due to any other unexplained reason. If not taken into account properly, this irregularity could result in the underestimation of brightness temperatures, when calibration information is smoothed over many scanlines. In *pygac*, this inconsistency is handled properly while calibration thermal channels.

In some cases it was found that, apart from the reset values, even the readings from any one of the four PRTs could also have very low suspicious values. This could also seriously affect the computation of brightness temperatures. *pygac* detects such anomalies and corrects them using interpolation of nearby valid PRT readings.




GAC I/O module
--------------

.. automodule:: pygac.gac_pod
   :members:
   :undoc-members:

The I/O module generates three HDF5 files, one containing reflectances, brightness temperatures, and lat/lon information. The other output file contains solar and satellite zenith and azimuth angles. And the third file contains quality flags.

The output file name format is:

``ECC_GAC_avhrr_satellitename_99999_yyyymmddThhmmsstZ_yyyymmddThhmmsstZ.h5``

and

``ECC_GAC_sunsatangles_satellitename_99999_yyyymmddThhmmsstZ_yyyymmddThhmmsstZ.h5``

and

``ECC_GAC_qualflags_satellitename_99999_yyyymmddThhmmsstZ_yyyymmddThhmmsstZ.h5``

where,

``ECC``: ESA CCI Clouds (This prefix can be changed/specified by the user)

``avhrr``: denoting that it contains reflectances and BTs

``sunsatangles``: denoting that it contains angles

``qualflags``: denoting that it contains quality flag information

``yyyymmddThhmmsstZ``: yy:year, mm:month, dd:day, hh:hour, mm:min, ss:sec, t:tenth of second (for the start and the end of the orbit).

Letters ``T`` and ``Z`` are separators for time info.

The value of ``99999`` is currently used instead of providing actual orbit number.

Appendices A, B and C provide detailed format of these files, including variable names, scaling, etc.


The start and end times in the header and in actual L1b data can be different for orbits. The mismatch can range from few milliseconds to days. It was decided to trust the time stamps in L1b data in *pygac*. After reorganizing based on scanline numbers (see issue highlighted above), the time stamps from the first and the last scanlines are taken as start and end times.


In some orbits, the latitude and longitude information contains corrupt values for a part/s of the orbit and these scanlines are not flagged in the corresponding scanline-by-scaline quality flags. Currently, *pygac* uses only a simple if_else construct to constrain valid range. Further improvement could be done using extra- and interpolation techniques.

For extremely warm and cold temperatures, the channel 3b is saturating producing irrelevant brightness temperatures. Such saturation is often not flagged in quality information. *pygac* currently uses a simple if_else construct to constrain valid range of Bts (170.0K<BT<350.0K). Such condition is also applied to split-window channels.


Supplement A: Structure of an output file containing reflectances and brightness temperatures
--------------

.. automodule:: pygac.gac_klm
   :members:
   :undoc-members:
   
Input: L1b file: ``NSS.GHRR.NN.D06279.S1800.E1955.B0711012.GC``

Output::

  Group how {
     variables:
       char channel_list(6, 9);
         :_lastModified = "2014-10-29T13:55:05Z";

     how:yaw_error = 0.0; // double
     how:roll_error = 0.0; // double
     how:pich_error = 0.0; // double
     how:startepochs = 1160157611L; // long
     how:endepochs = 1160164510L; // long
     how:platform = "noaa18";
     how:instrument = "avhrr";
     how:orbit_number = 99999; // int
     how:software = "pyGAC";
     how:version = "1.0";
   }

   Group image1 {
     variables:
       short data(13686, 409);
         :_lastModified = "2014-10-29T13:55:05Z";

     Group how {
       how:sun_earth_distance_correction_applied = "TRUE";
       how:sun_earth_distance_correction_factor = 0.9982412208987179; // double
     }

     Group what {
       what:product = "SATCH";
       what:quantity = "REFL";
       what:dataset_name = "Channel 1 reflectance";
       what:units = "%";
       what:gain = 0.01f; // float
       what:offset = 0.0f; // float
       what:missingdata = -32001; // int
       what:nodata = -32001; // int
       what:starttime = "180011";
       what:endtime = "195510";
       what:startdate = "20061006";
       what:enddate = "20061006";
     }

     image1:channel = "1";
     image1:description = "AVHRR ch1";
   }

   Group image2 {
     variables:
       short data(13686, 409);
         :_lastModified = "2014-10-29T13:55:05Z";

     Group how {
       how:sun_earth_distance_correction_applied = "TRUE";
       how:sun_earth_distance_correction_factor = 0.9982412208987179; // double
     }

     Group what {
       what:product = "SATCH";
       what:quantity = "REFL";
       what:dataset_name = "Channel 2 reflectance";
       what:units = "%";
       what:gain = 0.01f; // float
       what:offset = 0.0f; // float
       what:missingdata = -32001; // int
       what:nodata = -32001; // int
       what:starttime = "180011";
       what:endtime = "195510";
       what:startdate = "20061006";
       what:enddate = "20061006";
     }

     image2:channel = "2";
     image2:description = "AVHRR ch2";
   }

   Group image3 {
     variables:
       short data(13686, 409);
         :_lastModified = "2014-10-29T13:55:05Z";

     Group how {
     }

     Group what {
       what:product = "SATCH";
       what:quantity = "TB";
       what:dataset_name = "Channel 3b brightness temperature";
       what:units = "K";
       what:gain = 0.01f; // float
       what:offset = 273.15f; // float
       what:missingdata = -32001; // int
       what:nodata = -32001; // int
       what:starttime = "180011";
       what:endtime = "195510";
       what:startdate = "20061006";
       what:enddate = "20061006";
     }

     image3:description = "AVHRR ch3b";
     image3:channel = "3b";
   }

   Group image4 {
     variables:
       short data(13686, 409);
         :_lastModified = "2014-10-29T13:55:05Z";

     Group how {
     }

     Group what {
       what:product = "SATCH";
       what:quantity = "TB";
       what:dataset_name = "Channel 4 brightness temperature";
       what:units = "K";
       what:gain = 0.01f; // float
       what:offset = 273.15f; // float
       what:missingdata = -32001; // int
       what:nodata = -32001; // int
       what:starttime = "180011";
       what:endtime = "195510";
       what:startdate = "20061006";
       what:enddate = "20061006";
     }

     image4:channel = "4";
     image4:description = "AVHRR ch4";
   }

   Group image5 {
     variables:
       short data(13686, 409);
         :_lastModified = "2014-10-29T13:55:05Z";

     Group how {
     }

     Group what {
       what:product = "SATCH";
       what:quantity = "TB";
       what:dataset_name = "Channel 5 brightness temperature";
       what:units = "K";
       what:gain = 0.01f; // float
       what:offset = 273.15f; // float
       what:missingdata = -32001; // int
       what:nodata = -32001; // int
       what:starttime = "180011";
       what:endtime = "195510";
       what:startdate = "20061006";
       what:enddate = "20061006";
     }

     image5:channel = "5";
     image5:description = "AVHRR ch5";
   }

   Group image6 {
     variables:
       short data(13686, 409);
         :_lastModified = "2014-10-29T13:55:05Z";

     Group how {
       how:sun_earth_distance_correction_applied = "TRUE";
       how:sun_earth_distance_correction_factor = 0.9982412208987179; // double
     }

     Group what {
       what:product = "SATCH";
       what:quantity = "REFL";
       what:dataset_name = "Channel 3a reflectance";
       what:units = "%";
       what:gain = 0.01f; // float
       what:offset = 0.0f; // float
       what:missingdata = -32001; // int
       what:nodata = -32001; // int
       what:starttime = "180011";
       what:endtime = "195510";
       what:startdate = "20061006";
       what:enddate = "20061006";
     }

     image6:channel = "3a";
     image6:description = "AVHRR ch3a";
   }

   Group what {
     what:object = "SATP";
     what:sets = 6; // int
     what:version = "H5rad ?.?";
     what:date = "20061006";
     what:time = "180011";
   }

   Group where {

     Group lat {
       variables:
         int data(13686, 409);
           :_lastModified = "2014-10-29T13:55:05Z";

       Group what {
         what:dataset_name = "Latitude";
         what:units = "Deg";
         what:gain = 0.001f; // float
         what:offset = 0.0f; // float
         what:missingdata = -32001; // int
         what:nodata = -32001; // int
         what:starttime = "180011";
         what:endtime = "195510";
         what:startdate = "20061006";
         what:enddate = "20061006";
       }
     }

     Group lon {
       variables:
         int data(13686, 409);
           :_lastModified = "2014-10-29T13:55:05Z";

       Group what {
         what:dataset_name = "Longitude";
         what:units = "Deg";
         what:gain = 0.001f; // float
         what:offset = 0.0f; // float
         what:missingdata = -32001; // int
         what:nodata = -32001; // int
         what:starttime = "180011";
         what:endtime = "195510";
         what:startdate = "20061006";
         what:enddate = "20061006";
       }
     }

     where:num_of_pixels = 409; // int
     where:num_of_lines = 13686; // int
     where:xscale = 0.0f; // float
     where:yscale = 0.0f; // float
   }
  }



Supplement B: Structure of an output file containing Sun and satellite positions
--------------

.. automodule:: pygac.gac_klm
   :members:
   :undoc-members:

Input: L1b file: ``NSS.GHRR.NN.D06279.S1800.E1955.B0711012.GC``

Output::

  Group how {
     how:yaw_error = 0.0; // double
     how:roll_error = 0.0; // double
     how:pich_error = 0.0; // double
     how:startepochs = 1160157611L; // long
     how:endepochs = 1160164510L; // long
     how:platform = "noaa18";
     how:instrument = "avhrr";
     how:orbit_number = 99999; // int
     how:software = "pyGAC";
     how:version = "1.0";
   }

   Group image1 {
     variables:
       short data(13686, 409);
         :_lastModified = "2014-10-29T13:55:06Z";

     Group what {
       what:product = "SUNZ";
       what:quantity = "DEG";
       what:dataset_name = "Solar zenith angle";
       what:units = "Deg";
       what:gain = 0.01f; // float
       what:offset = 0.0f; // float
       what:missingdata = -32001; // int
       what:nodata = -32001; // int
       what:starttime = "180011";
       what:endtime = "195510";
       what:startdate = "20061006";
       what:enddate = "20061006";
     }

     image1:description = "Solar zenith angle";
   }

   Group image2 {
     variables:
       short data(13686, 409);
         :_lastModified = "2014-10-29T13:55:06Z";

     Group what {
       what:product = "SATZ";
       what:quantity = "DEG";
       what:dataset_name = "Satellite zenith angle";
       what:units = "Deg";
       what:gain = 0.01f; // float
       what:offset = 0.0f; // float
       what:missingdata = -32001; // int
       what:nodata = -32001; // int
       what:starttime = "180011";
       what:endtime = "195510";
       what:startdate = "20061006";
       what:enddate = "20061006";
     }

     image2:description = "Satellite zenith angle";
   }

   Group image3 {
     variables:
       short data(13686, 409);
         :_lastModified = "2014-10-29T13:55:06Z";

     Group what {
       what:product = "SSAZD";
       what:quantity = "DEG";
       what:dataset_name = "Relative satellite-sun azimuth angle";
       what:units = "Deg";
       what:gain = 0.01f; // float
       what:offset = 0.0f; // float
       what:missingdata = -32001; // int
       what:nodata = -32001; // int
       what:starttime = "180011";
       what:endtime = "195510";
       what:startdate = "20061006";
       what:enddate = "20061006";
     }

     image3:description = "Relative satellite-sun azimuth angle";
   }

   Group image4 {
     variables:
       short data(13686, 409);
         :_lastModified = "2014-10-29T13:55:06Z";

     Group what {
       what:product = "SUNA";
       what:quantity = "DEG";
       what:dataset_name = "Solar azimuth angle";
       what:units = "Deg";
       what:gain = 0.01f; // float
       what:offset = 180.0f; // float
       what:missingdata = -32001; // int
       what:nodata = -32001; // int
       what:starttime = "180011";
       what:endtime = "195510";
       what:startdate = "20061006";
       what:enddate = "20061006";
     }

     image4:description = "Solar azimuth angle";
   }

   Group image5 {
     variables:
       short data(13686, 409);
         :_lastModified = "2014-10-29T13:55:06Z";

     Group what {
       what:product = "SATA";
       what:quantity = "DEG";
       what:dataset_name = "Satellite azimuth angle";
       what:units = "Deg";
       what:gain = 0.01f; // float
       what:offset = 180.0f; // float
       what:missingdata = -32001; // int
       what:nodata = -32001; // int
       what:starttime = "180011";
       what:endtime = "195510";
       what:startdate = "20061006";
       what:enddate = "20061006";
     }

     image5:description = "Satellite azimuth angle";
   }

   Group what {
     what:object = "SATP";
     what:sets = 5; // int
     what:version = "H5rad ?.?";
     what:date = "20061006";
     what:time = "180011";
   }

   Group where {

     Group lat {
       variables:
         int data(13686, 409);
           :_lastModified = "2014-10-29T13:55:06Z";

       Group what {
         what:dataset_name = "Latitude";
         what:units = "Deg";
         what:gain = 0.001f; // float
         what:offset = 0.0f; // float
         what:missingdata = -32001; // int
         what:nodata = -32001; // int
         what:starttime = "180011";
         what:endtime = "195510";
         what:startdate = "20061006";
         what:enddate = "20061006";
       }
     }

     Group lon {
       variables:
         int data(13686, 409);
           :_lastModified = "2014-10-29T13:55:06Z";

       Group what {
         what:dataset_name = "Longitude";
         what:units = "Deg";
         what:gain = 0.001f; // float
         what:offset = 0.0f; // float
         what:missingdata = -32001; // int
         what:nodata = -32001; // int
         what:starttime = "180011";
         what:endtime = "195510";
         what:startdate = "20061006";
         what:enddate = "20061006";
       }
     }

     where:num_of_pixels = 409; // int
     where:num_of_lines = 13686; // int
     where:xscale = 0.0f; // float
     where:yscale = 0.0f; // float
   }
  }


Supplement C: Structure of an output file containing quality flags
--------------

.. automodule:: pygac.gac_klm
   :members:
   :undoc-members:

The file that contains quality flags has following information.

1) There is a variable called ``qual_flags/data``. It will have a dimension of (X,7), where X is the number of data records in the GAC orbit.

The 7 columns contain the following information.

=======  ==========================================================================================
Col 1    Scan line number
Col 2    Fatal error flag (scan line should not be used for analysis).
Col 3    Insufficient data for calibration (scan line should not be used for analysis).
Col 4    Insufficient data for navigation (scan line should not be used for analysis).
Col 5-7  whether solar contamination of blackbody occurred in in channels 3, 4, and 5 respectively.
=======  ==========================================================================================

If the values for these flags are greater than zero, then the data should not be used. If everything is normal, then all values should be zero.

2) There are also two important attributes that provide "last scan line number" and "total number of data records".

By combining information from column 1 and these two attributes,  the user is able to figure out where the gap occurs and also exact time for each scan line.
