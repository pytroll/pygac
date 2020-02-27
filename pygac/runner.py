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

"""Processing utilities for GAC/LAC KLM/POD files.

Functions:
    process_file: allows to process a given input file
    get_reader_class: allows to select the appropriate reader 
                      class for a given file.
"""

import datetime
import logging

from pygac.gac_klm import GACKLMReader
from pygac.gac_pod import GACPODReader
from pygac.lac_klm import LACKLMReader
from pygac.lac_pod import LACPODReader
from pygac.utils import file_opener


LOG = logging.getLogger(__name__)

_reader_classes = [GACKLMReader, LACKLMReader, GACPODReader, LACPODReader]


def get_reader_class(filename, fileobj=None):
    """Return the reader class that can read the GAC/LAC KLM/POD file.

        Args:
            filename (str): Path to GAC/LAC file
            fileobj: An open file object to read from. (optional)
    """
    found_reader = False
    for index, Reader in enumerate(_reader_classes):
        if Reader.can_read(filename, fileobj=fileobj):
            LOG.debug("%s can read the file." % Reader.__name__)
            found_reader = True
            break
    if not found_reader:
        raise ValueError('Unable to read the file "%s"' % filename)
    # Move the Reader in front of _reader_classes. Chance is high that the
    # next file is of the same kind.
    _reader_classes.insert(0, _reader_classes.pop(index))
    return Reader


def process_file(filename, start_line, end_line, fileobj=None):
    """Read, calibrate and navigate NOAA AVHRR GAC/LAC POD/KLM data.
    It creates three hdf5 files in the output location given by the pygac
    config file. The three files contain the avhrr data, quality flags,
    and sunsatangles.

       Args:
            filename (str): Path to GAC/LAC file
            start_line (int): First scanline to be processed (0-based)
            end_line (int): Last scanline to be processed (0-based),
                            set to 0 for the last available scanline
            fileobj: An open file object to read from. (optional)
    """
    tic = datetime.datetime.now()
    LOG.info("Process file: %s", str(filename))
    # Keep the file open while searching for the reader class and later
    # creation of the instance.
    with file_opener(fileobj or filename) as open_file:
        Reader = get_reader_class(filename, fileobj=open_file)
        reader = Reader.fromfile(filename, fileobj=open_file)
        reader.save(start_line, end_line)
    LOG.info("Processing took: %s", str(datetime.datetime.now() - tic))
