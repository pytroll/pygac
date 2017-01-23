
#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2012, 2014 Abhay Devasthale and Martin Raspaud

# Author(s):

#   Abhay Devasthale <abhay.devasthale@smhi.se>
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

"""Read a gac file.

"""

import argparse
import logging
from datetime import datetime

logger = logging.getLogger("")
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)


class MyFormatter(logging.Formatter):
    converter = datetime.fromtimestamp
    def formatTime(self, record, datefmt=None):
        ct = self.converter(record.created)
        if datefmt:
            s = ct.strftime(datefmt)
        else:
            t = ct.strftime("%Y-%m-%d %H:%M:%S")
            s = "%s.%03d" % (t, record.msecs)
        return s


formatter = MyFormatter('[ %(levelname)s %(name)s %(asctime)s] %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


def check_file_version(filename):
    try:
        with open(filename) as fdes:
            data = fdes.read(3)
    except IOError as err:
        raise IOError('Failed to read GAC file: {0}'.format(err))
    if data in ["CMS", "NSS", "UKM", "DSS"]:
        from pygac.gac_klm import main
        return main
    else:
        from pygac.gac_pod import main
        return main


def str2scanline(string):
    """Convert string to scanline.

    Make sure, the scanline is not negative.

    Args:
        string (str): String to be converted

    Returns:
        int: Scanline
    """
    integer = int(string)
    if integer < 0:
        raise argparse.ArgumentTypeError('Scanlines must be >= 0')
    return integer


def validate_args(args):
    if args.start_line > args.end_line:
        raise ValueError('Start Scanline > End Scanline')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Read, calibrate and navigate NOAA AVHRR GAC data')
    parser.add_argument('filename', type=str, help='GAC file to be processed')
    parser.add_argument('start_line', type=str2scanline,
                        help='First scanline to be processed (0-based)')
    parser.add_argument('end_line', type=str2scanline,
                        help='Last scanline to be processed (0-based, '
                             'set to 0 for the last available scanline)')
    args = parser.parse_args()
    validate_args(args)

    reader = check_file_version(args.filename)
    reader(args.filename, args.start_line, args.end_line)
