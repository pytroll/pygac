
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
    with open(filename) as fdes:
        data = fdes.read(3)
    if data in ["CMS", "NSS", "UKM", "DSS"]:
        from pygac.gac_klm import main
        return main
    else:
        from pygac.gac_pod import main
        return main



if __name__ == "__main__":
    import sys
    try:
        filename = sys.argv[1]
	start_line = sys.argv[2]
	end_line = sys.argv[3] 
    except IndexError:
	print "Usage: gac_run <filename> <start scan line number> <end scan line number>"
	sys.exit(1)

    reader = check_file_version(filename)
    try:
        reader(filename, start_line, end_line)
    except ValueError:
        print "Value error"

