Legacy Output
=============

GAC I/O module
--------------

The I/O module generates three HDF5 files, one containing reflectances,
brightness temperatures, and lat/lon information. The other output file
contains solar and satellite zenith and azimuth angles. And the third file
contains quality flags.

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

The value of ``99999`` is currently used instead of providing actual orbit
number.

Appendices A, B and C provide detailed format of these files, including
variable names, scaling, etc.


The start and end times in the header and in actual L1b data can be different
for orbits. The mismatch can range from few milliseconds to days. It was
decided to trust the time stamps in L1b data in Pygac. After reorganizing
based on scanline numbers (see issue highlighted above), the time stamps
from the first and the last scanlines are taken as start and end times.


In some orbits, the latitude and longitude information contains corrupt
values for a part/s of the orbit and these scanlines are not flagged in the
corresponding scanline-by-scaline quality flags. Currently, Pygac uses only
a simple if_else construct to constrain valid range. Further improvement could
be done using extra- and interpolation techniques.

For extremely warm and cold temperatures, the channel 3b is saturating
producing irrelevant brightness temperatures. Such saturation is often not
flagged in quality information. Pygac currently uses a simple if_else
construct to constrain valid range of Bts (170.0K<BT<350.0K). Such condition
is also applied to split-window channels.


.. automodule:: pygac.gac_io
   :members:
   :undoc-members:


Supplement A: Structure of an output file containing reflectances and brightness temperatures
---------------------------------------------------------------------------------------------

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
--------------------------------------------------------------------------------

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
------------------------------------------------------------------


The file that contains quality flags has following information.

1) There is a variable called ``qual_flags/data``. It will have a dimension of
(X,7), where X is the number of data records in the GAC orbit. The 7 columns
contain the following information.

=======  ==========================================================================================
Col 1    Scan line number
Col 2    Fatal error flag (scan line should not be used for analysis).
Col 3    Insufficient data for calibration (scan line should not be used for analysis).
Col 4    Insufficient data for navigation (scan line should not be used for analysis).
Col 5-7  whether solar contamination of blackbody occurred in in channels 3, 4, and 5 respectively.
=======  ==========================================================================================

If the values for these flags are greater than zero, then the data should not
be used. If everything is normal, then all values should be zero.

2) There are also two important attributes that provide "last scan line number"
and "total number of data records".

By combining information from column 1 and these two attributes, the user is
able to figure out where the gap occurs and also exact time for each scan line.
