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
logging.getLogger("pygac").addHandler(logging.NullHandler())
