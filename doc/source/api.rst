The :mod:`pygac` API
====================

pyGAC interface consists of a number of python modules. The following lines summarize the purpose of each python module.


gac_run.py      : main interface to all modules. Decides which version of AVHRR, and accordingly calls other modules.

gac_klm.py      : Reads L1b GAC data from KLM series of satellites (NOAA-15 and later) and does most of the computations.

gac_pod.py      : Reads L1b GAC data from POD series of satellites (NOAA-14 and before) and does most of the computations.

calibrate_klm.py: calibrates GAC data from KLM series, sets all necessary coefficients.

calibrate_pod.py: calibrates GAC data from POD series, sets all necessary coefficients.

gac_io.py       : I/O part of the software, writes two output files in HDF5 format.

geotiepoints.py : PyTROLL module used for interpolation.

astronomy.py    : PyTROLL module used for solar zenith.


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

ECC_avhrrGAC_noaaxx_yyyymmddZhhmmss_yyyymmddZhhmmss.h5
and
ECC_sunsatGAC_noaaxx_yyyymmddZhhmmss_yyyymmddZhhmmss.h5

where,

ECC: ESA CCI Clouds (This prefix can be specified by the user)

avhrrGAC: denoting that it contains reflectances and BTs

sunsatGAC: denoting that it contains angles

yyyymmddZhhmmss: year, month, day, hour, min, sec (for the start and the end of the orbit).

