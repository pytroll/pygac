Usage
-----

Through Satpy
~~~~~~~~~~~~~

The preferred way of using Pygac is through `Satpy`_. Results are returned as
dask-friendly `xarray`_ DataArrays with proper dataset/coordinate names and
additional metadata. It is also possible to select a user-defined range of
scanlines. Furthermore, Satpy provides many options for resampling,
visualizing and saving the data.

.. code-block:: python

    import satpy

    # Channel set for KLM satellites. For POD satellites the channels are
    # ['1', '2', '3', '4', '5'].
    channels = ['1', '2', '3a', '3b', '4', '5']
    ancillary = ['solar_zenith_angle',
                 'sensor_zenith_angle',
                 'solar_azimuth_angle',
                 'sensor_azimuth_angle',
                 'sun_sensor_azimuth_difference_angle',
                 'qual_flags',
                 'latitude',
                 'longitude']

    scene = satpy.Scene(filenames=['NSS.GHRR.NP.D15361.S0121.E0315.B3547172.SV'],
                        reader='avhrr_l1b_gaclac',
                        reader_kwargs={'tle_dir': '/path/to/tle/',
                                       'tle_name': 'TLE_%(satname)s.txt'})
    scene.load(channels + ancillary)


For a list of Satpy reader keyword arguments see `satpy.readers.avhrr_l1b_gaclac`_
and for further Pygac reader keyword arguments see :class:`pygac.reader.Reader`.
Especially it is possible to choose a different version of calibration
coefficients or even specify your own.

.. _Satpy: https://satpy.readthedocs.io
.. _xarray: https://xarray.pydata.org
.. _satpy.readers.avhrr_l1b_gaclac:
    https://satpy.readthedocs.io/en/stable/api/satpy.readers.avhrr_l1b_gaclac.html?highlight=avhrr_l1b_gaclac
.. _example notebook:
    https://github.com/pytroll/pytroll-examples/blob/main/satpy/avhrr_l1b_gaclac.ipynb


Direct Usage
~~~~~~~~~~~~

Alternatively you can also use Pygac directly.

.. code-block:: python

    from pygac import get_reader_cls
    
    filename = 'NSS.GHRR.NP.D15361.S0121.E0315.B3547172.SV'
    reader_cls = get_reader_cls(filename)
    reader = reader_cls(tle_dir='/path/to/tle', tle_name='TLE_%(satname)s.txt')
    reader.read(filename)

    channels = reader.get_calibrated_channels()
    lons, lats = reader.get_lonlat()
    scanline_times = reader.get_times()
    bad_quality_lines = reader.mask


Legacy CLI
~~~~~~~~~~

.. note::

    Usage of the legacy command line program ``pygac-run`` is deprecated in
    favour of the above options.

There is also a legacy command line program ``pygac-run`` which saves the
results to HDF5 and requires a configuration file.

Copy the template file ``etc/pygac.cfg.template`` to ``pygac.cfg`` and place
it in a directory as you please. Set the environment variable ``PYGAC_CONFIG_FILE``
pointing to the file. e.g.

.. code-block:: bash
 
  PYGAC_CONFIG_FILE=/home/user/pygac.cfg; export PYGAC_CONFIG_FILE

Also adapt the configuration file to your needs. The ``tledir`` parameter should
be set to where your Two Line Element (TLE) files are located.

Then call ``pygac-run`` on a GAC/LAC file.

.. code-block:: bash

  pygac-run testdata/NSS.GHRR.NL.D02187.S1904.E2058.B0921517.GC 0 0
 
The last two digits are the start and end scanline numbers, thus specifying the
portion of the GAC orbit that user wants to process. The first scanline number
starts at 0. If zeroes are specified at both locations, then the entire orbit
will be processed.

The result will be three hdf5 files, one with the calibrated AVHRR data,
the other with sun-satellite viewing geometry data and this third with
scanline quality information.


