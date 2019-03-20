#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2013, 2014 Martin Raspaud

# Author(s):

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

"""Test the GAC KLM reader.
"""

from pygac import gac_klm
import unittest
import subprocess as sp

# test_file = "test_data/NSS.GHRR.NL.D01183.S1313.E1458.B0400708.GC"
# ref_result = "test_data/ECC_avhrrGAC_noaa16_20010702Z131338_20010702Z145809.h5"
# my_result = "ECC_avhrrGAC_noaa16_20010702Z131338_20010702Z145809.h5"

# ref_sunsat = "test_data/ECC_sunsatGAC_noaa16_20010702Z131338_20010702Z145809.h5.orig"
# my_sunsat = "ECC_sunsatGAC_noaa16_20010702Z131338_20010702Z145809.h5"

test_file = "test_data/NSS.GHRR.NL.D02187.S1904.E2058.B0921517.GC"
ref_result = "test_data/ECC_GAC_sunsatangles_noaa16_99999_20020706T1904020Z_20020706T2058095Z.h5"
my_result = "/tmp/ECC_GAC_sunsatangles_noaa16_99999_20020706T1904020Z_20020706T2058095Z.h5"

ref_sunsat = "test_data/ECC_GAC_sunsatangles_noaa16_99999_20020706T1904020Z_20020706T2058095Z.h5"
my_sunsat = "/tmp/ECC_GAC_sunsatangles_noaa16_99999_20020706T1904020Z_20020706T2058095Z.h5"


class TestKLM(unittest.TestCase):

    def test_global(self):
        gac_klm.main(test_file, 0, 0)

        child = sp.Popen(["h5diff", ref_sunsat, my_sunsat],
                         stdout=sp.PIPE)
        streamdata = child.communicate()[0]
        retc = child.returncode
        self.assertTrue(retc == 0, msg=streamdata)

        child = sp.Popen(["h5diff", ref_result, my_result],
                         stdout=sp.PIPE)
        streamdata = child.communicate()[0]
        retc = child.returncode
        self.assertTrue(retc == 0, msg=streamdata)


if __name__ == '__main__':
    unittest.main()
