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

"""l1c files processor

"""
import datetime
import logging

from pygac.reader import ReaderError
from pygac.gac_klm import GACKLMReader
from pygac.gac_pod import GACPODReader
from pygac.lac_klm import LACKLMReader
from pygac.lac_pod import LACPODReader
from pygac.utils import file_opener


LOG = logging.getLogger(__name__)

_Readers = [GACKLMReader, LACKLMReader, GACPODReader, LACPODReader]

def get_reader(filename, fileobj=None):
    """Read the GAC/LAC KLM/POD data.

        Args:
            filename (str): Path to GAC/LAC file
            fileobj: An open file object to read from. (optional)
        """
    found_reader = False
    with file_opener(fileobj or filename) as open_file:
        for i, Reader in enumerate(_Readers):
            try:
                reader = Reader.fromfile(filename, fileobj=open_file)
                found_reader = True
                index = i
            except ReaderError as exception:
                LOG.debug("%s failed to read the file! %s" 
                          % (Reader.__name__, repr(exception)))
            finally:
                open_file.seek(0)
            if found_reader:
                break
    if not found_reader:
        raise ValueError('Unable to read the file "%s"' % filename)
    # Move the Reader in front of _Readers. Chance is high that the 
    # next file is of the same kind.
    _Readers.insert(0, _Readers.pop(index))
    return reader


def l1c_processor(filename, start_line, end_line, fileobj=None):
    """Level 1c file processor

       Args:
            filename (str): Path to GAC/LAC file
            start_line (int): First scanline to be processed (0-based)
            end_line (int): Last scanline to be processed (0-based),
                            set to 0 for the last available scanline
            fileobj: An open file object to read from. (optional)
    """
    tic = datetime.datetime.now()
    LOG.info("Process file: %s", str(filename))
    reader = get_reader(filename, fileobj=fileobj)
    reader.save(start_line, end_line)
    LOG.info("Processing took: %s", str(datetime.datetime.now() - tic))
