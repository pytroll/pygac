Usage
-----

Copy the template file ``etc/pygac.cfg.template`` to ``pygac.cfg`` and place
it in a directory as you please. Set the environment variable PYGAC_CONFIG_FILE
pointing to the file. E.g.::
 
  $> PYGAC_CONFIG_FILE=/home/user/pygac.cfg; export PYGAC_CONFIG_FILE

Also adapt the configuration file to your needs. The ``tledir`` parameter should
be set to where your Two Line Element (TLE) files are located.

The main script is ``pygac-run``. It automatically checks for the type of file
format and invokes either ``gac_pod.py`` (POD family, up to and including NOAA-14)
or ``gac_klm.py`` (KLM family, NOAA-15 and onwards including Metop-A and -B).

You can test it directly on the testdata included in the package. The result
will be three hdf5 files, one with the calibrated AVHRR data, the other with
sun-satellite viewing geometry data and this third with scanline quality
information::

 $> pygac-run testdata/NSS.GHRR.NL.D02187.S1904.E2058.B0921517.GC 0 0
 
The last two digits are the start and end scanline numbers, thus specifying the
portion of the GAC orbit that user wants to process.  The first scanline number
starts at 0. If zeroes are specified at both locations, then the entire orbit
will be processed.

