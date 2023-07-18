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
"""The setup module."""

from setuptools import setup
import os

if __name__ == '__main__':
    with open("README.md", "r") as fd:
        long_description = fd.read()

    requirements = ['docutils>=0.3',
                    'numpy>=1.8.0',
                    'pyorbital>=v0.3.2',
                    'h5py>=2.0.1',
                    'scipy>=0.8.0',
                    'python-geotiepoints>=1.1.8',
                    'bottleneck>=1.0.0']
    extras_require = {
        'dev': ['pytest', 'pre-commit', 'flake8']
    }

    setup(name='pygac',
          description='NOAA AVHRR GAC/LAC reader and calibration',
          author='The Pytroll Team',
          author_email='pytroll@googlegroups.com',
          classifiers=["Development Status :: 4 - Beta",
                       "Intended Audience :: Science/Research",
                       "License :: OSI Approved :: GNU General Public License v3 " +
                       "or later (GPLv3+)",
                       "Operating System :: OS Independent",
                       "Programming Language :: Python",
                       "Topic :: Scientific/Engineering"],
          url="https://github.com/pytroll/pygac",
          long_description=long_description,
          long_description_content_type='text/markdown',
          license='GPLv3',

          packages=['pygac'],
          package_data={'pygac': ['data/calibration.json']},

          # Project should use reStructuredText, so ensure that the docutils get
          # installed or upgraded on the target machine
          install_requires=requirements,
          extras_require=extras_require,
          scripts=[os.path.join('bin', item) for item in os.listdir('bin')],
          data_files=[('etc', ['etc/pygac.cfg.template']),
                      ('gapfilled_tles', ['gapfilled_tles/TLE_noaa16.txt'])],
          python_requires='>=3.8',
          zip_safe=False
          )
