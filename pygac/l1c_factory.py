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
    
    def __call__(filename, start_line, end_line, fileobj=None):
        tic = datetime.datetime.now()
        reader = self.Reader.fromfile(filename, fileobj=fileobj)
        self._save(reader)
        LOG.info("pygac took: %s", str(datetime.datetime.now() - tic))
    
    def _save(self, reader):
        channels = reader.get_calibrated_channels()
        sat_azi, sat_zen, sun_azi, sun_zen, rel_azi = reader.get_angles()
        qual_flags = reader.get_qual_flags()
        if (np.all(reader.mask)):
            LOG.error("All data is masked out. Stop processing!")
            raise ValueError("All data is masked out.")
        save_gac(
            reader.spacecraft_name, reader.utcs,
            reader.lats, reader.lons,
            channels[:, :, 0], channels[:, :, 1],
            channels[:, :, 2], channels[:, :, 3],
            channels[:, :, 4], channels[:, :, 5],
            sun_zen, sat_zen, sun_azi, sat_azi, rel_azi,
            qual_flags, start_line, end_line,
            reader.filename, reader.meta_data
        )

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
    def get_builder(cls, filename):
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
        # TODO: How do we determine if we need the lac_klm or lac_pod builder?
        if data in ["CMS", "NSS", "UKM", "DSS"]:
            return cls.__builders["gac_klm"]
        else:
            return cls.__builders["gac_pod"]
