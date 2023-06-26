## Version 1.7.2 (2023/06/26)


### Pull Requests Merged

#### Bugs fixed

* [PR 120](https://github.com/pytroll/pygac/pull/120) - Compatibility with numpy v1.24

#### Features added

* [PR 122](https://github.com/pytroll/pygac/pull/122) - Don't use deprecated distutils module.

#### Documentation changes

* [PR 123](https://github.com/pytroll/pygac/pull/123) - Update supported data formats in documentation ([2494](https://github.com/pytroll/satpy/issues/2494))

In this release 3 pull requests were closed.


## Version 1.7.1 (2022/12/09)


### Pull Requests Merged

#### Bugs fixed

* [PR 119](https://github.com/pytroll/pygac/pull/119) - Fix the documentation to include FRAC and light cleanup

#### Documentation changes

* [PR 119](https://github.com/pytroll/pygac/pull/119) - Fix the documentation to include FRAC and light cleanup

In this release 2 pull requests were closed.


## Version 1.7.0 (2022/10/31)

### Issues Closed

* [Issue 112](https://github.com/pytroll/pygac/issues/112) - Usage documentation uses wrong function name for 'get_reader_class' ([PR 113](https://github.com/pytroll/pygac/pull/113) by [@mraspaud](https://github.com/mraspaud))
* [Issue 78](https://github.com/pytroll/pygac/issues/78) - Handle completely missing ICT or space counts ([PR 117](https://github.com/pytroll/pygac/pull/117) by [@mraspaud](https://github.com/mraspaud))
* [Issue 22](https://github.com/pytroll/pygac/issues/22) - Thermal calibration error ([PR 117](https://github.com/pytroll/pygac/pull/117) by [@mraspaud](https://github.com/mraspaud))

In this release 3 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 117](https://github.com/pytroll/pygac/pull/117) - Fix calibration when 3b is totally deactivated ([78](https://github.com/pytroll/pygac/issues/78), [22](https://github.com/pytroll/pygac/issues/22))
* [PR 116](https://github.com/pytroll/pygac/pull/116) - Fix typos and deprecations
* [PR 113](https://github.com/pytroll/pygac/pull/113) - Fix usage typo ([112](https://github.com/pytroll/pygac/issues/112))

#### Features added

* [PR 118](https://github.com/pytroll/pygac/pull/118) - Bump up python versions
* [PR 115](https://github.com/pytroll/pygac/pull/115) - Add support for frac data

In this release 5 pull requests were closed.

## Version 1.6.0 (2022/08/04)

### Issues Closed

* [Issue 107](https://github.com/pytroll/pygac/issues/107) - Fix sun-earth distance correction
* [Issue 104](https://github.com/pytroll/pygac/issues/104) - Links to POD/KLM user guides are mixed up ([PR 105](https://github.com/pytroll/pygac/pull/105) by [@mraspaud](https://github.com/mraspaud))
* [Issue 103](https://github.com/pytroll/pygac/issues/103) - API is missing on readthedocs
* [Issue 40](https://github.com/pytroll/pygac/issues/40) - Calculated sun-earth distance correction factor once and add as attribute

In this release 4 issues were closed.

### Pull Requests Merged

#### Features added

* [PR 109](https://github.com/pytroll/pygac/pull/109) - Add new unpublished METOPC VIS calibration coefficients from Patmos-x

#### Documentation changes

* [PR 105](https://github.com/pytroll/pygac/pull/105) - Fix the pod and klm guide links ([104](https://github.com/pytroll/pygac/issues/104))

In this release 2 pull requests were closed.

## Version 1.5.0 (2022/01/10)

### Issues Closed

* [Issue 94](https://github.com/pytroll/pygac/issues/94) - Method get_tle_lines() raises index error when tle lines are empty strings
* [Issue 90](https://github.com/pytroll/pygac/issues/90) - Change in masked scanlines
* [Issue 87](https://github.com/pytroll/pygac/issues/87) - Unit test discovery ([PR 98](https://github.com/pytroll/pygac/pull/98) by [@sfinkens](https://github.com/sfinkens))
* [Issue 85](https://github.com/pytroll/pygac/issues/85) - Update documentation
* [Issue 80](https://github.com/pytroll/pygac/issues/80) - Reduce rounding error in POD reader adjust clock drift ([PR 84](https://github.com/pytroll/pygac/pull/84) by [@carloshorn](https://github.com/carloshorn))
* [Issue 76](https://github.com/pytroll/pygac/issues/76) - IndexError in PODReader._adjust_clock_drift ([PR 84](https://github.com/pytroll/pygac/pull/84) by [@carloshorn](https://github.com/carloshorn))
* [Issue 74](https://github.com/pytroll/pygac/issues/74) - Remove Cython dependency ([PR 83](https://github.com/pytroll/pygac/pull/83) by [@carloshorn](https://github.com/carloshorn))
* [Issue 46](https://github.com/pytroll/pygac/issues/46) - Pypi and presentation work
* [Issue 43](https://github.com/pytroll/pygac/issues/43) - Links to POD/KLM user guides are broken
* [Issue 39](https://github.com/pytroll/pygac/issues/39) - Update release in preparation for conda-forge package

In this release 10 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 91](https://github.com/pytroll/pygac/pull/91) - Invert bit position in quality flags
* [PR 84](https://github.com/pytroll/pygac/pull/84) - Refactor clock drift ([82](https://github.com/pytroll/pygac/issues/82), [81](https://github.com/pytroll/pygac/issues/81), [80](https://github.com/pytroll/pygac/issues/80), [76](https://github.com/pytroll/pygac/issues/76))
* [PR 79](https://github.com/pytroll/pygac/pull/79) - Fix tests on bigendian platforms
* [PR 75](https://github.com/pytroll/pygac/pull/75) - Typo

#### Features added

* [PR 98](https://github.com/pytroll/pygac/pull/98) - Switch to pytest for unit tests ([87](https://github.com/pytroll/pygac/issues/87))
* [PR 92](https://github.com/pytroll/pygac/pull/92) - Allow PathLike objects
* [PR 89](https://github.com/pytroll/pygac/pull/89) - Add github workflow
* [PR 86](https://github.com/pytroll/pygac/pull/86) - add METOP C coefficients
* [PR 84](https://github.com/pytroll/pygac/pull/84) - Refactor clock drift ([82](https://github.com/pytroll/pygac/issues/82), [81](https://github.com/pytroll/pygac/issues/81), [80](https://github.com/pytroll/pygac/issues/80), [76](https://github.com/pytroll/pygac/issues/76))
* [PR 83](https://github.com/pytroll/pygac/pull/83) - Replace cython filter ([74](https://github.com/pytroll/pygac/issues/74))

In this release 10 pull requests were closed.

## Version 1.4.0 (2020/08/06)

### Issues Closed

* [Issue 66](https://github.com/pytroll/pygac/issues/66) - Typos in calibration coefficients ([PR 67](https://github.com/pytroll/pygac/pull/67))
* [Issue 62](https://github.com/pytroll/pygac/issues/62) - Computation of Earth scene radiance ([PR 58](https://github.com/pytroll/pygac/pull/58))
* [Issue 60](https://github.com/pytroll/pygac/issues/60) - Improve readability of quality indicators bit unpacking ([PR 72](https://github.com/pytroll/pygac/pull/72))
* [Issue 57](https://github.com/pytroll/pygac/issues/57) - channel 4 BT to radiance conversion ([PR 67](https://github.com/pytroll/pygac/pull/67))
* [Issue 54](https://github.com/pytroll/pygac/issues/54) - Function check_file_version should not be part of pygac-run ([PR 55](https://github.com/pytroll/pygac/pull/55))
* [Issue 47](https://github.com/pytroll/pygac/issues/47) - Fix reading of renamed files

In this release 6 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 73](https://github.com/pytroll/pygac/pull/73) - Fix azimuth encoding
* [PR 67](https://github.com/pytroll/pygac/pull/67) - Correct coefficients ([66](https://github.com/pytroll/pygac/issues/66), [57](https://github.com/pytroll/pygac/issues/57))

#### Features added

* [PR 72](https://github.com/pytroll/pygac/pull/72) - Quality indicators ([60](https://github.com/pytroll/pygac/issues/60))
* [PR 71](https://github.com/pytroll/pygac/pull/71) - Expose new metadata
* [PR 70](https://github.com/pytroll/pygac/pull/70) - Remove config dependency from reader class
* [PR 58](https://github.com/pytroll/pygac/pull/58) - export coefficients to json file ([62](https://github.com/pytroll/pygac/issues/62))
* [PR 55](https://github.com/pytroll/pygac/pull/55) - Refactor gac-run ([54](https://github.com/pytroll/pygac/issues/54))

In this release 7 pull requests were closed.

## Version v1.3.1 (2020/02/07)

### Issues Closed

* [Issue 52](https://github.com/pytroll/pygac/issues/52) - Allow gzip input files ([PR 53](https://github.com/pytroll/pygac/pull/53))
* [Issue 49](https://github.com/pytroll/pygac/issues/49) - Calibration Coeffs Patmos-X 2017 ([PR 50](https://github.com/pytroll/pygac/pull/50))

In this release 2 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 50](https://github.com/pytroll/pygac/pull/50) - Fix typo in MetOp-B calibration ([49](https://github.com/pytroll/pygac/issues/49))
* [PR 48](https://github.com/pytroll/pygac/pull/48) - Update metadata when reading lonlat also

#### Features added

* [PR 53](https://github.com/pytroll/pygac/pull/53) - Allow gzip compressed files as input for klm and pod readers read method ([52](https://github.com/pytroll/pygac/issues/52))
* [PR 51](https://github.com/pytroll/pygac/pull/51) - Use the data format name to id ARS headered files

In this release 4 pull requests were closed.


## Version 1.3.0 (2019/12/05)


### Pull Requests Merged

#### Features added

* [PR 45](https://github.com/pytroll/pygac/pull/45) - Add LAC support ([5](https://github.com/pytroll/pygac/issues/5))
* [PR 42](https://github.com/pytroll/pygac/pull/42) - Update documentation
* [PR 41](https://github.com/pytroll/pygac/pull/41) - Add meta_data dictionary to reader class.

In this release 3 pull requests were closed.


## Version 1.2.1 (2019/11/15)

### Pull Requests Merged

#### Bugs fixed

* [PR 37](https://github.com/pytroll/pygac/pull/37) - Fixing geotiepoints attribute error for python2.7
* [PR 36](https://github.com/pytroll/pygac/pull/36) - Fix update of missing scanlines

#### Features added

* [PR 38](https://github.com/pytroll/pygac/pull/38) - Fix tests for python 3

In this release 3 pull requests were closed.

## Version 1.2.0 (2019/10/17)

### Issues Closed

* [Issue 33](https://github.com/pytroll/pygac/issues/33) - Make use of pytroll geotiepoints instead of using a deprecated copy
* [Issue 7](https://github.com/pytroll/pygac/issues/7) - Project URL points to wrong domain

In this release 2 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 30](https://github.com/pytroll/pygac/pull/30) - Feature update angles computation

#### Features added

* [PR 35](https://github.com/pytroll/pygac/pull/35) - Changed azimuth angle range to ]-180, 180]
* [PR 34](https://github.com/pytroll/pygac/pull/34) - Use the geotiepoints library
* [PR 32](https://github.com/pytroll/pygac/pull/32) - Updated documentation about azimuth angles
* [PR 31](https://github.com/pytroll/pygac/pull/31) - Refactor I/O

In this release 5 pull requests were closed.

## Version 1.1.0 (2019/06/12)

### Issues Closed

* [Issue 23](https://github.com/pytroll/pygac/issues/23) - Add support for Python3
* [Issue 20](https://github.com/pytroll/pygac/issues/20) - IndexError if orbit time beyond TLE range ([PR 21](https://github.com/pytroll/pygac/pull/21))
* [Issue 18](https://github.com/pytroll/pygac/issues/18) - Error in pygac data zenith angle
* [Issue 16](https://github.com/pytroll/pygac/issues/16) - Unit test failure

In this release 4 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 27](https://github.com/pytroll/pygac/pull/27) - Add get_times to angles computation
* [PR 24](https://github.com/pytroll/pygac/pull/24) - Fix unit tests
* [PR 21](https://github.com/pytroll/pygac/pull/21) - Fix TLE line identification ([20](https://github.com/pytroll/pygac/issues/20))
* [PR 17](https://github.com/pytroll/pygac/pull/17) - Python 3 compatibiity
* [PR 15](https://github.com/pytroll/pygac/pull/15) - Fix selection of user defined scanlines #2
* [PR 14](https://github.com/pytroll/pygac/pull/14) - Fix gac-run
* [PR 9](https://github.com/pytroll/pygac/pull/9) - Fixes: Dependency and docs
* [PR 3](https://github.com/pytroll/pygac/pull/3) - Feature clock fixed nan tsm scanline timestamp

#### Features added

* [PR 28](https://github.com/pytroll/pygac/pull/28) - Simplify TLE computation
* [PR 26](https://github.com/pytroll/pygac/pull/26) - Python-3 Compatibility
* [PR 25](https://github.com/pytroll/pygac/pull/25) - Fix style
* [PR 19](https://github.com/pytroll/pygac/pull/19) - Add support for TIROS-N
* [PR 13](https://github.com/pytroll/pygac/pull/13) - Update installation.rst with TLE files
* [PR 11](https://github.com/pytroll/pygac/pull/11) - Add scanline timestamps to qualflags file
* [PR 10](https://github.com/pytroll/pygac/pull/10) - Usability improvements
* [PR 6](https://github.com/pytroll/pygac/pull/6) - Add new attributes to /how
* [PR 4](https://github.com/pytroll/pygac/pull/4) - Improve interface of gac_run.py

In this release 17 pull requests were closed.
