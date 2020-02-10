#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014, 2019 Pygac Developers

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

"""Factory for l1c files.

"""
import logging

from pygac.utils import file_opener

LOG = logging.getLogger(__name__)

class L1cFactory(object):
    """Factory for level 1c files."""
    def __call__(self, filename, start_line, end_line):
        # check file version
        processor = self.check_file_version
        # process
        processor(filename, start_line, end_line)
    
    @staticmethod
    def check_file_version(filename):
        try:
            with file_opener(filename) as fdes:
                try:
                    data = fdes.read(3).decode()
                except UnicodeDecodeError:
                    data = None
        except IOError as err:
            raise IOError('Failed to read GAC file: {0}'.format(err))
        # TODO: discuss where this information belongs to.
        # Maybe every reader needs a can_read class method
        # which checks if this file belongs to the reader
        # Then the factory loops over the reader classes to determine
        # who can read it (would require to preload the classes)
        # which seems to be avoided by the developers
        # At least, this desission does not belong into the run script!
        if data in ["CMS", "NSS", "UKM", "DSS"]:
            from pygac.gac_klm import main
            return main
        else:
            from pygac.gac_pod import main
            return main
