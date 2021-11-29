Methods
=======

Calibration
-----------

At present, calibration coefficients provided by Andrew Heidinger
(NOAA) under SCOPE-CM project are applied for all satellites. The current
version is *PATMOS-x, v2017r1*, including provisional coefficients for MetOp-C.

The solar channel calibration (Channels 1 and 2, and Channel 3a if available)
takes into account inter-satellite differences and is derived using
amalgamation of different calibration references including the most recent
MODIS Collection 6 data, in-situ targets, and simultaneous nadir
observations. The detailed methodology is presented in Heidinger et al.
(2010). The resulting (inter)calibration coefficients are of highest quality
and follow the Global Climate Observing System (GCOS) standards for
deriving fundamental climate data records.

The reflectances are normalized by a factor (as a function of day of a year)
to account for changing Earth-Sun distance. However, it is left to the
user to apply further normalization using cosine of solar zenith
angle (depending on application in question).

The thermal channel intercalibration is done from scratch, starting from
obtaining Platinum Resistance Thermometer (PRT), space and Internal
Calibration Target (ICT, blackbody) counts. ICT temperatures are obtained
from PRT counts based on coefficients provided in the POD and KLM Data User
Guides (Kidwell, 2000). For each thermal channel, a smoothing window of 51
successive PRT, ICT and space counts is used to obtain robust gain values and
to dampen undue high frequency fluctuations in the count data (Trishchenko,
2002). This window size is easily configurable to suit user needs. Section
7.1.2.5 of `KLM User Guide`_ presents the summary of equations implemented
in Pygac to calibrate thermal channels, including non-linearity correction
(Walton et al. 1998).

.. _KLM User Guide:
    https://www.ncei.noaa.gov/pub/data/satellite/publications/podguides/TIROS-N%20thru%20N-14/

The PRT readings are supposed to be present in a specific order, for example,
the reset value followed by four readings from the PRTs. However, this may
not always be the case for orbits that contain data gaps or due to any other
unexplained reason. If not taken into account properly, this irregularity
could result in the underestimation of brightness temperatures, when
calibration information is smoothed over many scanlines. In Pygac, this
inconsistency is handled properly while calibrating thermal channels.

In some cases it was found that, apart from the reset values, even the
readings from any one of the four PRTs could also have very low suspicious
values. This could also seriously affect the computation of brightness
temperatures. Pygac detects such anomalies and corrects them using
interpolation of nearby valid PRT readings.


Geolocation
-----------

Each GAC row has total 409 pixels. But lat-lon values are provided for every
eigth pixel starting from pixel 5 and ending at pixel 405. Using Numpy, Scipy
and Pyresample packages, inter- and extrapolation is carried out to obtain
geolocation for each pixel along the scanline.

If the GAC data belongs to POD family, then clock drift errors are used to
adjust existing Lat-Lon information. Here, Pygac makes use of `PyOrbital`_
package. Pygac interpolates the clock offset and adjusts the nominal scan
times to the actual scan times. Since the geolocation was computed using the
nominal scan times, Pygacinterpolates the latitudes and longitudes to the
actual scan times using spherical linear interpolation, aka slerp. However,
in the case of a clock drift error greater than the scan rate of the dataset,
the latitude and longitude for each pixel of the scan lines that cannot have
an interpolated geolocation (typically at the start or end of the dataset)
are recomputed. This is done using pyorbital, which in turn uses TLEs to
compute the position of the satellite at each scan time and the instrument
geometry compute the longitude and latitude of each pixel of the dataset.
Since this operation can be quite costly, the interpolation is preferred
whenever possible.

.. _PyOrbital:
    https://pyorbital.readthedocs.io


Computation of Angles
---------------------

The azimuth angles are calculated using `get_alt_az`_ and `get_observer_look`_
from pyorbital. The azimuth described in the link is measured as clockwise
from North instead of counter-clockwise from South. Counter clockwise from
south would be the standard for a right-handed orthogonal coordinate system.
Pygac was updated to use the same definition for angles as pyorbital (2019,
September, version > 1.1.0). Previous versions used azimuth +/-180 degrees,
which correspond to degrees clockwise from south. All angles are converted to
degrees. All azimuth angles are converted to range ]-180, 180] (2019 October
version > 1.1.0 ). Note that ]-180, 180] is an open interval.


.. _get_alt_az:
    https://pyorbital.readthedocs.io/en/latest/#pyorbital.astronomy.get_alt_az
.. _get_observer_look:
    https://pyorbital.readthedocs.io/en/latest/#pyorbital.orbital.Orbital.get_observer_look


Correction of Satellite Location
--------------------------------

Whenever possible, Pygac uses RPY corrections along with other orbital
parameters to compute accurate satellite location (e.g. instead of assuming
constant altitude). However, RPY corrections are not available for all NOAA
satellites. In case of the majority of the POD family satellites, these
corrections are set to zero.


Correction of Scanline Timestamps
---------------------------------

The geolocation in Pygac depends on accurate scanline timestamps. However,
these may be corrupt, especially for older sensors. Assuming a constant
scanning rate, Pygac attempts to fix them using extrapolation based on the scan
line number and a reference time.

Finding the right reference time is difficult due to the multitude of
possible timestamp corruptions. But the combination of the following three
options proved to be a robust reference in many situations:
Timestamp of the first scanline, median time offset of all scanlines and header
timestamp. See
:meth:`pygac.reader.Reader.correct_times_median` and
:meth:`pygac.reader.Reader.correct_times_thresh`
for details.

Finally, not only timestamps but also scanline numbers may be corrupt.
Therefor lines with erroneous scanline numbers are removed before
extrapolation, see :meth:`pygac.reader.Reader.correct_scan_line_numbers`.


Scan-Motor-Issue
----------------

Between 2001 and 2004 GAC data from NOAA-14, NOAA-15, and NOAA-16 frequently
contain a significant amount of noise towards an edge of the swath. As
reported by `Schlundt et al (2017)`_, section 5.2, this is probably caused by a
temporary scan-motor issue. Pygac tries to identify and mask affected pixels.

.. _Schlundt et al (2017):
    https://climate.esa.int/media/documents/Cloud_Technical-Report-AVHRR-GAC-FCDR-generation_v1.0.pdf
