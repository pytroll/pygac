Changelog
=========

%%version%% (unreleased)
------------------------

- Update changelog. [Martin Raspaud]

- Bump version: 1.0.0 → 1.0.1. [Martin Raspaud]

- Use version number from version.py. [Martin Raspaud]

- Update changelog. [Martin Raspaud]

- Bump version: 0.1.0 → 1.0.0. [Martin Raspaud]

- Do style cleanup for gac_pod. [Martin Raspaud]

- Allow PYGAC_CONFIG_FILE to be missing. [Martin Raspaud]

- Merge branch 'master' into develop. [Adam.Dybbroe]

- Merge branch 'pre-master' [Adam.Dybbroe]

- Merge branch 'feature-clock' into pre-master. [Adam.Dybbroe]

- Time stamp for mid-morning orbit investigated. [Abhay Devasthale]

- Time stamp for mid-night orbit investigated. [Abhay Devasthale]

- The support for NOAA-6 and -8 added. [Abhay Devasthale]

- Calibration support for NOAA-6 and -8 added. [Abhay Devasthale]

- Merge branch 'feature-clock' into pre-master. [Martin Raspaud]

- Changed noaa07 and noaa09 to noaa7 and noaa9 resp. [Abhay Devasthale]

- I/O part updated. PyGAC will first cut the start and end of the orbit
  if invalid geolocation is present and THEN apply user provided start
  and end scanline numbers to output the orbit. [Abhay Devasthale]

- Absolute value of minimum scanline number used to remove negative scan
  line numbers. [Abhay Devasthale]

- Modifications to cut out orbits that contain invalid lat/lon values at
  the start and the end of the tracks. [Abhay Devasthale]

- NOAA-11 IndexError messages fixed. Scan line number format changed
  from unsigned to signed integer. [Abhay Devasthale]

- Update to handle out of range values of time stamps. [Abhay
  Devasthale]

- Updated launch dates to match Andy's. Also changes calibration
  coefficients for MetOps. [Abhay Devasthale]

- Bug fix while subsetting channel 2. [Abhay Devasthale]

- Better handling of incorrect time stamps in the L1b data. Correct
  reorganizing of scanlines if the orbit is missing starting scanlines.
  Minor bug fixes in the header reading. [Abhay Devasthale]

- Merge branch 'feature-clock' of github.com:adybbroe/pygac into
  feature-clock. [Martin Raspaud]

- Changed handeling of buggy lat,lons and scanline numbers (Abhay)
  [Nina.Hakansson]

- Bugfix: masking invalid brightness temperatures. [Nina.Hakansson]

  Before values set to -32001 where multiplied with 100 and saved as in16
  and got actual values -32768.


- Merge branch 'feature-clock' of ssh://github.com/adybbroe/pygac into
  feature-clock. [Nina.Hakansson]

- Rpy fix for pod-satellites (Abhay) [Nina.Hakansson]

- Fixed the unit tests. [Martin Raspaud]

- Don't use pod's fixed_error_correction rpy parameters. [Martin
  Raspaud]

  Use klm's constant_[roll,pitch,yaw]_attitude_error parameters instead.

- Fix the rpy problem for older pod scenes. [Martin Raspaud]

- Pep8ify. [Adam Dybbroe]

- Updated lauch dates for several satellites. [Sara.Hornquist]

- Quick fix to handle channel 3a/3b switch. [Nina.Hakansson]

  Maybe better to not interpolate.


- Corrected launch date for metopb. [Sara.Hornquist]

- Channel 3a/3b transition set to missing data. [Nina.Hakansson]

- Adding metopb to gac_klm.py. [Nina.Hakansson]

- Bugfix: apply earth-sun distance correction. [Nina.Hakansson]

- If all pixels will be masked out, stop processing; i.e. do not produce
  output files. [Sara.Hornquist]

- When using startline/endline, set times in the filenames that
  corresponds to the actual start and end time of the selected part of
  the orbit. Do not use the line numbers in the filenames any more.
  (changes from Abhay) [Sara.Hornquist]

- Merge branch 'feature-clock' of ssh://github.com/adybbroe/pygac into
  feature-clock. [Sara.Hornquist]

- Remove obsolete hdf5 files. [Adam Dybbroe]

- Bug fix, one attribute was misspelled. [Sara.Hornquist]

- Corrections in writing products, in order to match PPS v2014: - Make
  sure the string attribute has got the right format. - Changed nodata
  used for lat/lon. (the old nodata was inside valid range)
  [Sara.Hornquist]

- Fix readme file. [Adam Dybbroe]

- Gain factor applied to lat lon values. [Abhay Devasthale]

- Minor correction in calibration file. NOAA-7 amd MetOp-A tested.
  [Abhay Devasthale]

- Channel 3 BTs over the Antarctica corrected. [Abhay Devasthale]

- Minor changes to output file name. [Abhay Devasthale]

- Filenaming and quality flag changes. [Abhay Devasthale]

  * Filenaming convention is changed to harmonize with pps convention.
  * Quality flag file updated.


- Merge branch 'feature-clock' of github.com:adybbroe/pygac into
  feature-clock. [Abhay Devasthale]

  Conflicts:
  	pygac/gac_pod.py


- Fix the clock drift correction. [Martin Raspaud]

- Remove testing 2.6 in travis until scipy dependency is removed.
  [Martin Raspaud]

- Try another hdf5 dev package. [Martin Raspaud]

- Fix to allow h5py to build on 2.6. [Martin Raspaud]

- Don't raise an error if PYGAC_CONFIG_FILE doesn't point to a file.
  [Martin Raspaud]

- Fixing travis for 2.6. [Martin Raspaud]

- Recent changes to GAC IO. [Abhay Devasthale]

- Bugfixing. [Martin Raspaud]

- Bugfix calibration coefficients. [Martin Raspaud]

- Added missing calibration coefficients. [Martin Raspaud]

- Add the gac reader generic class. [Martin Raspaud]

- CI on 2.6 and add the PYGAC env var. [Martin Raspaud]

- Completing calibration coefficients. [Martin Raspaud]

- Finished factorizing, hopefully. [Martin Raspaud]

- Add slerp tests. [Martin Raspaud]

- Numpy 1.8.0 needed at least. [Martin Raspaud]

- Revamped tests. [Martin Raspaud]

- Implemented clock drift for pod. [Martin Raspaud]

- Add slerp computations. [Martin Raspaud]

- Add a simple clock drift adjustment (line shifting) [Martin Raspaud]

- WIP: Update calibration coeffs. [Martin Raspaud]

- Finish factorizing code for calibration. Some calibration coeffs
  missing. [Martin Raspaud]

- WIP: Clock drift and refactoring. [Martin Raspaud]

- Cleaning, and beginning of refactoring. [Martin Raspaud]

- Supplements A, B and C added. [abhaydd]

- Updating documentation. [abhaydd]

- Updating pygac api documentation. [abhaydd]

- Updated text on command-line usage. [abhaydd]

- Update usage.rst. [abhaydd]

- Update usage.rst. [abhaydd]

- Bugfix. [Adam Dybbroe]

- Added for scipy dependency. [Adam Dybbroe]

- Added requirements file, for Travis... [Adam Dybbroe]

- Added support for travis. [Adam Dybbroe]

- Added buttons on readme page for code health etc. [Adam Dybbroe]

- Added customization support for Landscape. [Adam Dybbroe]

- Smoothing window for thermal channel calibration adjusted.
  [Abhay.Devasthale]

- Updates on time information in output files. No 10th seconds, and
  seconds-since-1970 is now properly set. [Sara.Hornquist]

- Merge branch 'pre-master' of github.com:adybbroe/pygac into pre-
  master. [Sara.Hornquist]

- Dumping of debugging info on screen is avoided in gac_pod.py.
  [Abhay.Devasthale]

- Update in output files: attribute what/time do not have tenth-of-
  second any more. [Sara.Hornquist]

- Updated documentation on filenames. [Sara.Hornquist]

- Negative reflectances replaced by MISSING_DATA. [Abhay.Devasthale]

- Replaced nighttime reflectances with MISSING_DATA. [Abhay.Devasthale]

- POD: Refined the tle search to get the nearest match. [Martin Raspaud]

  In the case of old satellites, the tle data can be quite scarse. For that
  reason, the find_tle_index function was enhanced to provide the closest
  match to the required date.

- Bugfix in pod, and cleanup. [Martin Raspaud]

  - A correct determination of which sensor was generating each prt has been
    implemented, allowing the data to miss scanlines. It is based on the
    scanline numbers provided in the data
  - The pod data is also cleaned up before after reading.
  - The code has been cleaned up a little, to follow python standards.

- Remove astronomy.py, depend on pyorbital instead. [Martin Raspaud]

- Added h5py as a requirement in setup. [Adam Dybbroe]

- Merge branch 'pre-master' of github.com:adybbroe/pygac into pre-
  master. [Adam Dybbroe]

- Add some test scripts, and remove test data. [Martin Raspaud]

- Added documentation. [Abhay.Devasthale]

- Update api.rst. [abhaydd]

- Update api.rst. [abhaydd]

- Update api.rst. [abhaydd]

- Update api.rst. [abhaydd]

- Minor editorial. [Adam Dybbroe]

- Fixing Manifest and setup. [Adam Dybbroe]

- Updated usage docs. [Adam Dybbroe]

- Adding a bit of documentation and the test case. [Adam Dybbroe]

- Add empty (sphinx) docs. [Adam Dybbroe]

- Adding configuration and logging. [Adam Dybbroe]

- Merge branch 'master' into develop. [Adam Dybbroe]

- Changed readme. [Adam Dybbroe]

- Making a python package out of it. [Adam Dybbroe]

- Initial commit. [Adam Dybbroe]


