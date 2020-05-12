#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2020 Pygac Developers

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

"""Module configuration class.

Read and manage module configuration.
"""
import logging
import os
import sys
try:
    import configparser
except ImportError:
    import ConfigParser as configparser

if sys.version_info.major < 3:
    class FileNotFoundError(OSError):
        pass

LOG = logging.getLogger(__name__)


class Configuration(configparser.ConfigParser, object):
    """Configuration container for pygac."""

    config_file = ''

    def read(self, config_file):
        """Read and parse the configuration file

        Args:
            config_file (str): path to the config file

        Note:
            In contrast to the parent method, this implementation
            only accepts a single file and raises an exception if
            the given file is not readable.
        """
        if not os.path.isfile(config_file):
            raise FileNotFoundError(
                'Given config path "%s" is not a file!' % config_file)
        try:
            super(Configuration, self).read(config_file)
        except configparser.NoSectionError:
            LOG.error('Failed reading configuration file: "%s"'
                      % config_file)
            raise
        self.config_file = config_file

    def get(self, *args, **kwargs):
        """python 2 compatibility for fallback attribute"""
        if sys.version_info.major < 3:
            if 'fallback' in kwargs:
                fallback = kwargs.pop('fallback')
            else:
                fallback = None
            try:
                value = super(Configuration, self).get(*args, **kwargs)
            except (configparser.NoSectionError, configparser.NoOptionError):
                value = fallback
        else:
            value = super().get(*args, **kwargs)
        return value


_config = Configuration()


def get_config(initialized=True):
    """Retrun the module configuration.

    Args:
        initialized (bool): if true ensure that the configuration has
                            been initialized (default: True)
    """
    global _config
    if initialized and not _config.sections():
        LOG.info('Configuration has not been initialized. Use'
                 ' environment variable "PYGAC_CONFIG_FILE"')
        try:
            config_file = os.environ["PYGAC_CONFIG_FILE"]
        except KeyError:
            LOG.error('Environment variable PYGAC_CONFIG_FILE not set!')
            raise
        _config.read(config_file)
    return _config


def read_config_file(config_file):
    """Read a given config file."""
    config = get_config(initialized=False)
    config.read(config_file)
