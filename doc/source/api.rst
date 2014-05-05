The :mod:`pygac` API
====================

Overview of the code


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

