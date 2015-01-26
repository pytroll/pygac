The :mod:`pygac` API
====================

PyGAC interface consists of a number of python modules. The following schematic shows a general structure and the processing flow of PyGAC.   

First it is determined whether the GAC data comes from the second (i.e. NOAA-14 and before) or the third generation (NOAA-15 and onwards) AVHRR instrument by :mod:`gac_run.py`.







GAC KLM reader
--------------

.. automodule:: pygac.gac_klm
   :members:
   :undoc-members:



GAC POD reader
--------------

.. automodule:: pygac.gac_pod
   :members:
   :undoc-members:


GAC I/O module
--------------

.. automodule:: pygac.gac_pod
   :members:
   :undoc-members:

The I/O module generates two HDF5 files, one contaning reflectances, brightness temperatures, and lat/lon information.
The other output file contains solar and satellite zenith and azimuth angles.

The output file name format is:

ECC_GAC_avhrr_satellitename_99999_yyyymmddThhmmsstZ_yyyymmddThhmmsstZ.h5

and

ECC_GAC_sunsatangles_satellitename_99999_yyyymmddThhmmsstZ_yyyymmddThhmmsstZ.h5

where,

ECC: ESA CCI Clouds (This prefix can be specified by the user)

avhrr: denoting that it contains reflectances and BTs

sunsatangles: denoting that it contains angles

yyyymmddThhmmsstZ: yy:year, mm:month, dd:day, hh:hour, mm:min, ss:sec, t:tenth of second (for the start and the end of the orbit).

Letters T and Z are separators for time info.

The value of 99999 is currently used instead of providing actual orbit number.
