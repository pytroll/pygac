Introduction
============

Supported Data Format
---------------------

Pygac reads AVHRR GAC (Global Area Coverage), LAC (Local Area Coverage) and FRAC
(Full Resolution Area Coverage) level 1b data from NOAA, which is described in
the `POD`_ (NOAA-14 and before) and `KLM`_ (NOAA-15 and following) user guides.
The data can be obtained from `NOAA CLASS`_, where you can also find a
comprehensive `introduction`_.

.. note::

    Pygac can only read AVHRR data with 10 bits/pixel. This can be specified in
    your NOAA CLASS order or user preferences (see screenshot below).

.. image:: images/noaa_class_preferences.png
    :width: 350
    :alt: Screenshot of AVHRR data extraction preferences

.. _NOAA CLASS:
    https://www.class.noaa.gov/
.. _POD:
    https://github.com/user-attachments/files/17407676/POD-Guide.pdf
.. _KLM:
    https://github.com/user-attachments/files/17361330/NOAA_KLM_Users_Guide.pdf
.. _introduction:
    https://www.class.noaa.gov/release/data_available/avhrr/index.htm


Supported Sensors
-----------------
Pygac currently supports AVHRR generations 1-3 onboard NOAA (TIROS-N, NOAA-6
and onwards) and MetOp satellites.


.. _here:
    https://www.avl.class.noaa.gov/release/data_available/avhrr/index.htm


Related Projects
----------------

- `pygac-fdr`_: Generate a fundamental data record of AVHRR GAC data using
  Pygac.
- `level1c4pps`_: Prepare AVHRR GAC data for `NWCSAF/PPS`_ using Pygac.

.. _level1c4pps: https://github.com/foua-pps/level1c4pps
.. _NWCSAF/PPS: https://www.nwcsaf.org/16
.. _pygac-fdr: https://github.com/pytroll/pygac-fdr
