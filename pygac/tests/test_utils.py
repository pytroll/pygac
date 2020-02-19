#tests/test_utils.py
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

import unittest
import io
import gzip

from pygac.utils import is_file_object, file_opener


class TestUtils(unittest.TestCase):
    """Test pygac.utils functions"""

    longMessage = True

    def setUp(self):
        """Set up the test."""
        self.normal_message = b'normal message'
        self.gzip_message_decoded = b'gzip message'
        with io.BytesIO() as f:
            with gzip.open(f, mode='wb') as g:
                g.write(self.gzip_message_decoded)
            f.seek(0)
            self.gzip_message_encoded = f.read()

    def test_is_file_object(self):
        """Test is_file_object function."""
        # true test
        with io.BytesIO(self.normal_message) as fileobj:
            self.assertTrue(is_file_object(fileobj))
        # false test
        self.assertFalse(is_file_object("test.txt"))
        # duck type test
        class Duck(object):
            def read(self, n):
                return n*b'\00'
            def seekable(self):
                return True
            def close(self):
                pass
        duck = Duck()
        self.assertTrue(is_file_object(duck))

    def test_file_opener(self):
        """Test file_openter function"""
        # On normal file (check also if it remains open)
        with io.BytesIO(self.normal_message) as f:
            with file_opener(f) as g:
                message = g.read()
            self.assertFalse(f.closed)
        self.assertEqual(message, self.normal_message)
        # On gzip file
        with io.BytesIO(self.gzip_message_encoded) as f:
            with file_opener(f) as g:
                message = g.read()
        self.assertEqual(message, self.gzip_message_decoded)


def suite():
    """Test suite for test_utils."""
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestUtils))
    return mysuite


if __name__ == '__main__':
    unittest.main()
