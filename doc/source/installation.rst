Installation
------------

pygac reads NOAA AVHRR Global Area Coverage (GAC) data, and performs state of
the art calibration and navigation.

You can install the latest stable release of the software via the python package index (pypi)::

  $> pip install --prefix=$PREFIX pygac

where ``$PREFIX`` denotes your desired installation prefix. If you want access to the full 
source code and/or test the latest version under development you can download the pygac 
source code from github_::

  $> git clone git://github.com/pytroll/pygac.git
  $> cd pygac

and then run::

  $> pip install --prefix=$PREFIX .

or, if you want to hack the package::

  $> pip install -e --prefix=$PREFIX .

To uninstall the package run::

  $> pip uninstall pygac


At runtime
----------

Let $PREFIX be the prefix of your pygac installation. To make the python 
interpreter aware of the pygac module, update the ``$PYTHONPATH`` environment
variable::

  $> export PYTHONPATH=$PREFIX/lib/python2.7/site-packages:$PYTHONPATH

Furthermore, update the ``$PATH`` environment variable to have the pygac 
scripts available in your session::

  $> export PATH=$PREFIX/bin:$PATH


.. _github: http://github.com/pytroll/pygac

TLE files
---------
The pygac package requires Two-Line Element files stored per-satellite
in files with names such as TLE_noaa19.txt. The name format and directory can be
configured in the config file (see the Usage section). The contents should be the
historical TLEs, i.e. a concatenation of just lines 1 and 2 without the satellite
name. These can be downloaded from space-track.org. Once logged in use URL such as
(replace CCCCC with the five-digit NORAD catalogue ID):

https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/1990-01-01--2029-01-01/NORAD_CAT_ID/CCCCC/orderby/EPOCH ASC/format/tle
