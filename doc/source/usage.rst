Usage
-----

Copy the template file *pygac.cfg_template* to *pygac.cfg* and place
it in a directory as you please. Set the environment variable PYGAC_CONFIG_FILE
pointing to the file. E.g.::
 
  $> PYGAC_CONFIG_FILE=/home/user/pygac.cfg; export PYGAC_CONFIG_FILE

Also adapt the configuration file to your needs. The *tledir* parameter should
be set to where your Two Line Element (TLE) files are located.

The main script is *gac_run.py*. It automatically checks for the type of file
file format and invokes either gac_pod.py (pre-KLM type satellites) or
gac_klm.py (KLMNN' satellites). You can test it directly on the testdata
included in the package. The result will be two hdf5 files, one with the
calibrated AVHRR data, and one with sun-satellite viewing geometry data::

 $> python pygac/gac_run.py testdata/NSS.GHRR.NL.D02187.S1904.E2058.B0921517.GC


