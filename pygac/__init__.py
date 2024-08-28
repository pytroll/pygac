#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c14526.ad.smhi.se>
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

try:
    from pygac.version import version as __version__  # noqa
except ModuleNotFoundError:
    raise ModuleNotFoundError(
        "No module named pygac.version. This could mean "
        "you didn't install 'pygac' properly. Try reinstalling ('pip "
        "install').")

import logging
from pygac.configuration import get_config, read_config_file  # noqa
from pygac.runner import get_reader_class, process_file  # noqa

# add a NullHandler to prevent messages in sys.stderr if the using application does
# not use logging, but pygac makes logging calls of severity WARNING and greater.
# See https://docs.python.org/3/howto/logging.html (Configuring Logging for a Library)
logging.getLogger('pygac').addHandler(logging.NullHandler())
