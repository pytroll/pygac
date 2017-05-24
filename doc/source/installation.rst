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
