Usage
-----

Copy the template file *pygac.cfg_template* to *pygac.cfg* and place
it in a directory as you please. Set the environment variable PYGAC_CONFIG_FILE
pointing to the file. E.g.::
 
  $> PYGAC_CONFIG_FILE=/home/user/pygac.cfg; export PYGAC_CONFIG_FILE

Also adapt the configuration file to your needs. The *tledir* parameter should
be set to where your Two Line Element (TLE) files are located.


A simple use case::

  >>> from pyspectral.rsr_read import RelativeSpectralResponse
  >>> from pyspectral.solar import (SolarIrradianceSpectrum, TOTAL_IRRADIANCE_SPECTRUM_2000ASTM)
  >>> modis = RelativeSpectralResponse('eos', '2', 'modis')
  >>> modis.load(channel='20', scale=0.001)
  >>> solar_irr = SolarIrradianceSpectrum(TOTAL_IRRADIANCE_SPECTRUM_2000ASTM, dlambda=0.005)
  >>> sflux = solar_irr.inband_solarflux(modis.rsr)
  >>> print("Solar flux over Band: ", sflux)
  ('Solar flux over Band: ', 2.002927764514423)

And, here is how to derive the solar reflectance (removing the thermal part) of
the Aqua MODIS 3.7 micron band::

  >>> from pyspectral.nir_reflectance import Calculator
  >>> sunz = 80.
  >>> tb3 = 290.0
  >>> tb4 = 282.0
  >>> refl37 = Calculator(modis.rsr, solar_flux=sflux)
  >>> print refl37.reflectance_from_tbs(sunz, tb3, tb4)
  0.251177702956
