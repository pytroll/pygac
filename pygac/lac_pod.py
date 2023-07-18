#!/usr/bin/python
# Copyright (c) 2014-2019.
#

# Author(s):

#   Abhay Devasthale <abhay.devasthale@smhi.se>
#   Adam Dybbroe <adam.dybbroe@smhi.se>
#   Sajid Pareeth <sajid.pareeth@fmach.it>
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

"""Reader for LAC POD data."""

import logging

import numpy as np

from pygac.pod_reader import PODReader, main_pod
from pygac.lac_reader import LACReader

LOG = logging.getLogger(__name__)

scanline = np.dtype([("scan_line_number", ">i2"),
                     ("time_code", ">u2", (3, )),
                     ("quality_indicators", ">u4"),
                     ("calibration_coefficients", ">i4", (10, )),
                     ("number_of_meaningful_zenith_angles_and_earth_location_appended",
                      ">u1"),
                     ("solar_zenith_angles", "i1", (51, )),
                     ("earth_location", ">i2", (102, )),
                     ("telemetry", ">u4", (35, )),
                     ("sensor_data", ">u4", (3414, )),
                     ("add_on_zenith", ">u2", (10, )),
                     ("clock_drift_delta", ">u2"),  # only in newest version
                     ("spare3", "u2", (337, ))])


class LACPODReader(LACReader, PODReader):
    """The LAC POD reader.

    The `scan_points` attributes provides the position of the longitude and latitude points to
    compute relative to the full swath width.

    The offset attribute tells where in the file the scanline data starts.
    """

    def __init__(self, *args, **kwargs):
        """Init the LAC POD reader."""
        LACReader.__init__(self, *args, **kwargs)
        self.scanline_type = scanline
        self.offset = 14800
        self.scan_points = np.arange(2048)


def main(filename, start_line, end_line):
    """Generate a l1c file."""
    return main_pod(LACPODReader, filename, start_line, end_line)


if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2], sys.argv[3])
