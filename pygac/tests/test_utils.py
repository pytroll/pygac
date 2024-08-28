#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013, 2014 Martin Raspaud

# Author(s):

#   Carlos Horn <carlos.horn@external.eumetsat.int>

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

"""Test pygac.utils module
"""

import gzip
import io
import os
import unittest
from unittest import mock

import numpy as np

from pygac.utils import file_opener, calculate_sun_earth_distance_correction


class TestUtils(unittest.TestCase):
    """Test pygac.utils functions"""

    longMessage = True

    @mock.patch('pygac.utils.open', mock.MagicMock(return_value=io.BytesIO(b'file content')))
    def test_file_opener_1(self):
        """Test if a file is redirected correctly through file_opener."""
        with file_opener('path/to/file') as f:
            content = f.read()
        self.assertEqual(content, b'file content')

    def test_file_opener_2(self):
        """Test file_opener with file objects and compression"""
        # prepare test
        normal_message = b'normal message'
        gzip_message_decoded = b'gzip message'
        with io.BytesIO() as f:
            with gzip.open(f, mode='wb') as g:
                g.write(gzip_message_decoded)
            f.seek(0)
            gzip_message_encoded = f.read()
        # on normal file (check also if it remains open)
        with io.BytesIO(normal_message) as f:
            with file_opener(f) as g:
                message = g.read()
            self.assertFalse(f.closed)
        self.assertEqual(message, normal_message)
        # on gzip file
        with io.BytesIO(gzip_message_encoded) as f:
            with file_opener(f) as g:
                message = g.read()
        self.assertEqual(message, gzip_message_decoded)

    @mock.patch('pygac.utils.open', mock.MagicMock(side_effect=FileNotFoundError))
    def test_file_opener_3(self):
        """Test file_opener with PathLike object"""
        # prepare test
        class RawBytes(os.PathLike):
            def __init__(self, filename, raw_bytes):
                self.filename = str(filename)
                self.raw_bytes = raw_bytes

            def __fspath__(self):
                return self.filename

            def open(self):
                return io.BytesIO(self.raw_bytes)

        filename = '/path/to/file'
        file_bytes = b'TestTestTest'
        test_pathlike = RawBytes(filename, file_bytes)
        with file_opener(test_pathlike) as f:
            content = f.read()
        self.assertEqual(content, file_bytes)

        # test with lazy loading open method (open only in context)
        class RawBytesLazy(RawBytes):
            def open(self):
                self.lazy_opener_mock = mock.MagicMock()
                self.lazy_opener_mock.__enter__.return_value = io.BytesIO(self.raw_bytes)
                return self.lazy_opener_mock

        test_pathlike = RawBytesLazy(filename, file_bytes)
        with file_opener(test_pathlike) as f:
            content = f.read()
        self.assertEqual(content, file_bytes)
        test_pathlike.lazy_opener_mock.__exit__.assert_called_once_with(None, None, None)

    def test_calculate_sun_earth_distance_correction(self):
        """Test function for the sun distance corretction."""
        corr = calculate_sun_earth_distance_correction(3)
        np.testing.assert_almost_equal(corr, 0.96660494, decimal=7)
