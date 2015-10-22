#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 Abhay Devasthale, Adam.Dybbroe, Martin Raspaud

# Author(s):

#   Adam.Dybbroe <a000680@c14526.ad.smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


try:
    with open("./README", "r") as fd:
        long_description = fd.read()
except IOError:
    long_description = ""


from setuptools import setup
import imp

version = imp.load_source('pygac.version', 'pygac/version.py')

setup(name='pygac',
      version="v0.1.0",
      description='NOAA AVHRR GAC reader and calibration',
      author='Abhay Devasthale',
      author_email='adam.dybbroe@smhi.se',
      classifiers=["Development Status :: 4 - Beta",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: GNU General Public License v3 " +
                   "or later (GPLv3+)",
                   "Operating System :: OS Independent",
                   "Programming Language :: Python",
                   "Topic :: Scientific/Engineering"],
      url="https://github.com/adybbroe/pygac",
      long_description=long_description,
      license='GPLv3',

      packages=['pygac'],

      # Project should use reStructuredText, so ensure that the docutils get
      # installed or upgraded on the target machine
      install_requires=['docutils>=0.3',
                        'numpy>=1.8.0', 'pyorbital>=v0.3.2',
                        'h5py'],
      extras_require={'geolocation interpolation': ['python-geotiepoints'],
                      },
      scripts=[],
      data_files=[('etc', ['etc/pygac.cfg.template']),
                  ('gapfilled_tles', ['gapfilled_tles/TLE_noaa16.txt'])],
      test_suite="pygac.tests.suite",
      tests_require=[],

      zip_safe=False
      )
