#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c14526.ad.smhi.se>

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

import logging
import os

from pygac.version import __version__  # noqa

LOG = logging.getLogger(__name__)
try:
    CONFIG_FILE = os.environ['PYGAC_CONFIG_FILE']
except KeyError:
    LOG.error('Environment variable PYGAC_CONFIG_FILE not set!')
    CONFIG_FILE = ''

if not os.path.exists(CONFIG_FILE) or not os.path.isfile(CONFIG_FILE):
    LOG.warning(str(CONFIG_FILE) + " pointed to by the environment "
                + "variable PYGAC_CONFIG_FILE is not a file or does not exist!")
