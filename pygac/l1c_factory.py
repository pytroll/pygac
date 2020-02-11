#!/usr/bin/env python

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
import datetime
import logging

from pygac.gac_klm import GACKLMReader
from pygac.gac_pod import GACPODReader
from pygac.lac_klm import LACKLMReader
from pygac.lac_pod import LACPODReader
from pygac.utils import file_opener
from pygac.gac_io import save_gac


LOG = logging.getLogger(__name__)

class L1cBuilder(object):
    def __init__(self, Reader):
        self.Reader = Reader
    
    def __call__(self, filename, start_line, end_line, fileobj=None):
        tic = datetime.datetime.now()
        LOG.info("Building file: %s", str(filename))
        reader = self.Reader.fromfile(filename, fileobj=fileobj)
        reader.save(start_line, end_line)
        LOG.info("Processing took: %s", str(datetime.datetime.now() - tic))

class L1cFactory(object):
    """Factory for level 1c files."""
    __builders = {
        "gac_klm": L1cBuilder(GACKLMReader),
        "gac_pod": L1cBuilder(GACPODReader),
        "lac_klm": L1cBuilder(LACKLMReader),
        "lac_pod": L1cBuilder(LACPODReader)
    }
    def __call__(self, filename, start_line, end_line, fileobj=None):
        # check which builder to use
        builder = self.get_builder(filename, fileobj=fileobj)
        # process
        builder(filename, start_line, end_line, fileobj=None)
    
    @classmethod
    def get_builder(cls, filename, fileobj=None):
        try:
            if fileobj is None:
                open_file = file_opener(filename)
            else:
                open_file = file_opener(fileobj)
            with open_file as fdes:
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
        # who can read it.
        # At least, this desission does not belong into the run script!
        # How do we determine if we need the lac_klm or lac_pod builder?
        if data in ["CMS", "NSS", "UKM", "DSS"]:
            builder = cls.__builders["gac_klm"]
        else:
            builder = cls.__builders["gac_pod"]
        return builder
        
