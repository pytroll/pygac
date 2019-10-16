#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author(s):

#   Stephan Finkensieper <stephan.finkensieper@dwd.de>

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

"""Test utilities."""

import numpy as np


class CalledWithArray(object):
    """Adapter for arrays in mock.assert_called_with()."""

    def __init__(self, array):
        self.array = array

    def __repr__(self):
        return repr(self.array)

    def __eq__(self, other):
        return np.all(self.array == other)
