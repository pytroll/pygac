"""Module configuration class.

Read and manage module configuration.
"""
import configparser
import logging
import os

LOG = logging.getLogger(__name__)


class Configuration(configparser.ConfigParser):
    """Configuration container for pygac."""

    config_file = ""

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
            super().read(config_file)
        except configparser.NoSectionError:
            LOG.error('Failed reading configuration file: "%s"'
                      % config_file)
            raise
        self.config_file = config_file


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
            LOG.error("Environment variable PYGAC_CONFIG_FILE not set!")
            raise
        _config.read(config_file)
    return _config


def read_config_file(config_file):
    """Read a given config file."""
    config = get_config(initialized=False)
    config.read(config_file)
