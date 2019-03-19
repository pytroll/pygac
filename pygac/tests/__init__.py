#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013, 2014 Martin Raspaud

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>

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

"""The tests package.
"""

from pygac.tests import test_calibrate_pod, test_slerp, test_calibrate_klm, \
    test_pod, test_corrections, test_reader, test_io
import unittest


def suite():
    """The global test suite.
    """
    mysuite = unittest.TestSuite()
    tests = (test_slerp, test_calibrate_klm, test_calibrate_pod,
             test_pod, test_corrections, test_reader, test_io)
    for test in tests:
        mysuite.addTests(test.suite())

    return mysuite
